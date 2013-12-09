(ns iiit
  (:import 
    [java.io Writer File OutputStream FileOutputStream BufferedReader BufferedWriter FileReader FileWriter]
    [java.util Collections])
  (:use 
    [clojure stacktrace]
    [clojure.set :only (union intersection difference)]
    [clojure.contrib repl-utils trace pprint str-utils (io :only (spit append-spit make-parents delete-file))]
    [clojure.contrib.shell :only (sh)]))

(println "Starting to load name-space iiit ..")

;(set! *warn-on-reflection* true)

(defrecord Residue [^String name ^String type ^String chain ^int id])

; The RNA graph is represented as (Adjacency list rep - persistent) :
; 1. nodes       : Map of Node names => Nodes
; 2. cbond-edges : Map of Node names => Tuple of nodes (prev cbond, next cbond)
; 3. hbond-edges : Map of Node names => Set of Node names

(defrecord RNAAsGraph [nodes cbond-edges hbond-edges])

(defn res-id      [RAG r] (:id ((:nodes RAG) r)))
(defn cbond-edges [RAG r] ((:cbond-edges RAG) r))
(defn hbond-edges [RAG r] ((:hbond-edges RAG) r))

(defn convert-atoms-pdb-to-res-seq
  "Forms a sequence of residues from the Atoms PDB file.
   The information taken out is: res-type, res-num & chain
   Each residue has a unique name: R_<chain><res_num>
   Additionally each residue has a unique id.
   Returns a sequence of Residue records (nodes) having the following keys:
   name (R_<chain><res_num>), type (A,U,G,C), chain, id"
  [^String atom-pdb-file-name]
  (with-open [rdr (BufferedReader. (FileReader. atom-pdb-file-name))]
    (let [residues
          (filter (comp not nil?) 
            (distinct
              (map (fn [l]
                     (when (> (count l) 30)
                       (let [res-type (re-gsub #"\s+" "" (subs l 17 20))
                             res-name (re-gsub #"\s+" "" (subs l 21 26))
                             chain    (subs res-name 0 1)]
                         {:name (str "R_" res-name) :type res-type :chain chain}))) 
                (line-seq rdr))))]
      (loop [residue-seq []
             res-seq (seq residues)
             id 0]
        (if res-seq
          (let [r (first res-seq)
                R (Residue. (:name r) (:type r) (:chain r) id)]
            (recur (conj residue-seq R) (next res-seq) (inc id)))
          residue-seq)))))

;(count (convert-atoms-pdb-to-res-seq "C:\\Docs\\ramsagar\\Input"))
;1505
;(count (convert-atoms-pdb-to-res-seq "C:\\Docs\\ramsagar\\1VQO-min.coor"))
;10549

(defn form-residue-graph-with-cbonds
  "Gets list of residues from the atom-pdb-file by calling convert-atoms-pdb-to-res-seq.
   Returns the - 
   1. Mapping from Node names => Nodes
   2. Cbond-Edges as : Map of Node names => Tuple of nodes [prev_cbond, next_cbond]
   Residue with id k, is convalently bonded to residue with id k-1, k+1
   unless the residue is at the beginning/end of the chain."
  [^String atom-pdb-file-name]
  (letfn [(f [residues prev-res nodes cbond-edges]
            (if-let [curr-res (first residues)]
              (let [curr-id (:name curr-res)
                    new-nodes (assoc nodes curr-id curr-res)
                    new-cbond-edges
                    (let [prev-id (:name prev-res)
                          chains-match? (= (:chain prev-res) (:chain curr-res))
                          new-cbond-edges
                          (assoc cbond-edges prev-id (conj (cbond-edges prev-id) (if chains-match? curr-id nil)))
                          new-cbond-edges
                          (assoc new-cbond-edges curr-id [(if chains-match? prev-id nil)])]
                      new-cbond-edges)]
                (recur (next residues) curr-res new-nodes new-cbond-edges))
              (let [new-cbond-edges
                    (assoc cbond-edges (:name prev-res) (conj (cbond-edges (:name prev-res)) nil))]
                (list nodes new-cbond-edges))))]
    (let [residues (convert-atoms-pdb-to-res-seq atom-pdb-file-name)]
      (if-let [r (first residues)]
        (f (next residues) r (assoc {} (:name r) r) (assoc {} (:name r) [nil]))
        (list nil nil)))))

(defn form-hbond-edge-map
  "Takes the hbonds file produced by the program Hbonds.java as input.
   Takes following info from the file - 
   1. Donor residue num, chain
   2. Acceptor residue num, chain
   *Note: do not change the format of the file as returned by Hbonds.java*
   Returns - 
   1. Hbond-Edges as : Map of Node names => Set of Node names 
  "
  [hbonds-file]
  (letfn [(put-hbonds [hbond-edges seq-strs]
            (if seq-strs
              (let [line (first seq-strs)
                    tokens (re-split #"\s+" line)]
                (if (= 1 (count(nth tokens 3)))
                  (let [first-res-id (str "R_" (nth tokens 5) (nth tokens 4) )
                        second-res-id (str "R_" (nth tokens 10) (nth tokens 9))
                        new-hbond-edges 
                        (assoc hbond-edges first-res-id
                          (if (contains? hbond-edges first-res-id)
                            (conj (hbond-edges first-res-id) second-res-id)
                            #{second-res-id}))
                        new-hbond-edges 
                        (assoc new-hbond-edges second-res-id
                          (if (contains? new-hbond-edges second-res-id)
                            (conj (new-hbond-edges second-res-id) first-res-id)
                            #{first-res-id}))]
                    (recur new-hbond-edges (next seq-strs)))
                  (recur hbond-edges (next seq-strs))))
              hbond-edges))]
    (put-hbonds {} (filter #(not= % "") (re-split #"\n" (slurp hbonds-file))))))

(defn form-RNA-graph
  "Forms the RNA graph taking an options map with the following 2 keys -
   1. atoms-pdb-file-name
   2. hbonds-file-name (Generated by Hbonds.java)

   The RNA graph is represented as (Adjacency list rep. - persistent) :
     1. nodes       : Map of Node names => Nodes
     2. cbond-edges : Map of Node names => Tuple of nodes (prev cbond, next cbond)
     3. hbond-edges : Map of Node names => Set of Node names"
  [opts]
  (let [{:keys [atoms-pdb-file-name hbonds-file-name]} opts
        [nodes cbond-edges] (form-residue-graph-with-cbonds atoms-pdb-file-name)
        hbond-edges         (form-hbond-edge-map hbonds-file-name)]
    (RNAAsGraph. nodes cbond-edges hbond-edges)))

;;;;;;;;;;;;;;;;;;;;;;; UTILITY ;;;;;;;;;;;;;;;;;;;;;;;; 
(defn get-hbonds-seq
  "Returns a sequence of all hbonds in the RAG.
   The single argument version returns all hbonds (undirected)
   The two argument version returns distinct directed hbonds"
  ([RAG]
    (mapcat (fn [r] (map #(vector r %) (hbond-edges RAG r))) (keys (:hbond-edges RAG))))
  ([RAG one-dir?]
    (filter (fn [[r1 r2]] (< (res-id RAG r1) (res-id RAG r2))) (get-hbonds-seq RAG))))

; (count (get-hbonds-seq RAG))
; (count (get-hbonds-seq RAG true))

(defn get-cbonds-seq
  "Returns a sequence of all cbonds in the RAG.
   The single argument version returns all cbonds (undirected)
   The two argument version returns distinct directed cbonds"
  ([RAG ]
    (filter (comp not nil?) 
      (mapcat (fn [r] (map #(when (not (nil? %)) (vector r %)) (cbond-edges RAG r))) (keys (:cbond-edges RAG)))))
  ([RAG one-dir?]
    (filter (fn [[r1 r2]] (< (res-id RAG r1) (res-id RAG r2))) (get-cbonds-seq RAG))))

; (count (get-cbonds-seq RAG))
; (count (get-cbonds-seq RAG true))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



(defn get-RNA-graph-stats 
  "Returns a map containing the following keys -
   - Number of Residue nodes
   - Number of Directed CBond edges
   - Number of Directed HBond edges"
  [RAG]
  {"Number of Residue nodes"        (count (:nodes RAG))
   "Number of Directed CBond edges" (count (get-cbonds-seq RAG true))
   "Number of Directed HBond edges" (count (get-hbonds-seq RAG true))})
  
(defn print-RNA-graph-to-file 
  "Prints a RAG in a file named by fileName.  The file is overwritten every time.
   It prints in order :
    * Residue Nodes (1 per line)
    * CBond edges (1 residue's per line)
    * Hbond edges (1 residue's per line)
  "
  [RAG ^String fileName]
  (with-open [w (BufferedWriter. (FileWriter. fileName))]
    (.write w "===================== Residue Nodes =====================\n")
    (doseq [k (sort (:nodes RAG))]    (.write w (str k "\n")))
    (.write w "\n===================== CBond Edges =====================\n")
    (doseq [k (sort (:cbond-edges RAG))] (.write w (str k "\n")))
    (.write w "\n===================== HBond Edges =====================\n")
    (doseq [k (sort (:hbond-edges RAG))] (.write w (str k "\n")))))


;;;;;;;;;;;;;;;;;;;;;;; N-Lets ;;;;;;;;;;;;;;;;;;;;;;;; 
(declare nlets-non-distinct nlets-distinct)

(defn nlets-non-distinct
  [RAG n]
  (if (<= n 1) 
    nil
    (if (= n 2)
      (get-hbonds-seq RAG)
      (mapcat 
        (fn [n-1let] 
          (mapcat 
            (fn [r] 
              (filter (comp not nil?) 
                (map 
                  (fn [h] 
                    (when-not (some #{h} n-1let) (conj n-1let h))) 
                  (hbond-edges RAG r)))) 
            n-1let))
        (nlets-distinct RAG (dec n))))))


(defn nlets-distinct
  "Returns all nlets of length 'n'
   Definition of nlet:
   A set of N connected (by hbonds) nodes.
   Each nlet is just a list of nodes
  "
  [RAG n]
  (set (map sort (nlets-non-distinct RAG n))))

(def nlets-distinct (memoize nlets-distinct))

; (dotimes [i 34] (println (+ 2 i) "," (count (nlets-distinct RAG (+ 2 i)))))
; (doseq [l (take 10 (nlets-distinct RAG 3))] (println l))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;; STEMS ;;;;;;;;;;;;;;;;;;;;;;;;;; 
(defn stem 
  "Forms a stem if possible of len:max-len starting from the first hbond in vec-s, vec-t"
  [RAG vec-s vec-t dir-s dir-t curr-len max-len]
  (let [last-s (vec-s (dec (count vec-s)))
        last-t (vec-t (dec (count vec-t)))
        next-s ((cbond-edges RAG last-s) dir-s)
        next-t ((cbond-edges RAG last-t) dir-t)]
    (when (and next-s next-t)
      (let [is-hbond? (contains? (hbond-edges RAG next-s) next-t)]
        (when is-hbond?
          (let [max-s (if (= 1 dir-s) next-s (vec-s 0))
                min-t (if (= 1 dir-t) (vec-t 0) next-t)]
            (when (< (res-id RAG max-s) (res-id RAG min-t))
              (let [new-len (inc curr-len)
                    new-vec-s (conj vec-s next-s)
                    new-vec-t (conj vec-t next-t)]
                (if (= max-len new-len)
                  [new-vec-s new-vec-t]
                  (recur RAG new-vec-s new-vec-t dir-s dir-t new-len max-len))))))))))

(defn stems
  "Get all stems of length len.

   Definition of a stem of length n:
   A set S of two sets S1, S2 of equal number of 'n'nodes such that, 
   S1 = {r1, r2, ... } with ri covalently bonded to both ri-1 and ri+1 ; 
   S2 = {t1, t2, ... } with ti covalently bonded to both ti-1 and ti+1 ; 
   additionally each ri and ti form hydrogen bonds.

   Returns the sequence of stems. Each stem is a tuple. First elem is one side, second element the other side. 
  "  
  [RAG len]
  (let [hbonds (get-hbonds-seq RAG true)]
    (filter (comp not nil?) 
      (mapcat (fn [h] (list (stem RAG [(h 0)] [(h 1)] 0 1 1 len))) hbonds))))

(def stems (memoize stems))

; (doseq [s (stems RAG 13)] (println s))
; (dotimes [i 25] (println (+ 2 i) "," (count (stems RAG (+ 2 i)))))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;; STEM LOOPS ;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-next-residues-from
  "Gets next residues in the chain starting from from-res and looks ahead num-res"
  [RAG from-res num-res]
  (loop [i num-res, s-loop [], prev-r from-res]
    (if (> i 0)
      (when-let [next-r ((cbond-edges RAG prev-r) 1)]
        (recur (dec i) (conj s-loop next-r) next-r))
      s-loop)))

(defn is-sloop-devoid-of-hbonds-among-itself?
  "A complicated condition for the loop in a stem to be devoid of hbonds"
  [RAG max-s s-loop min-t]
  (let [S-loop (concat [max-s] s-loop [min-t])
        intra-h-bonds-absent? 
        (nil? 
          (first
            (for [r1 S-loop, r2 S-loop 
                  :when (and (< (res-id RAG r1) (res-id RAG r2))
                          (not= r2 ((cbond-edges RAG r1) 1))
                          (contains? (hbond-edges RAG r1) r2)
                          (not (and (= r1 max-s) (= r2 min-t))))]
              [r1 r2])))]
    (if (= 2 (count s-loop))
      (and intra-h-bonds-absent? (not (contains? (hbond-edges RAG (s-loop 0)) (s-loop 1))))
      intra-h-bonds-absent?)))

(defn connect
  "Try to connect low residue with high residue with num-res nodes"
  [RAG low high num-res]
  (let [connection (get-next-residues-from RAG low num-res)
        is-connected? (and (not (nil? (seq connection)))
                        (= ((cbond-edges RAG high) 0) (connection (dec (count connection)))))]
    (when is-connected?
      connection)))

(defn stem-is-stem-loop
  "Checks if a stem is a stem loop with num-in-loop nodes in the loop"
  [RAG stem num-in-loop]
  (let [max-s ((stem 0) 0)
        min-t ((stem 1) 0)
        s-loop (connect RAG max-s min-t num-in-loop)
        is-s-loop-ok? (and s-loop (is-sloop-devoid-of-hbonds-among-itself? RAG max-s s-loop min-t))]
    (when is-s-loop-ok?
      s-loop)))

(defn stem-loops ;aka hair-pin loops
  "Get stemLoops of length lenStem and numInLoop residues in the loop. A.k.a hairPinLoops

   Definition of a stemLoop :
   Stem loops are stems of length lenStem with an additional property that r1 and t1 are on the same chain 
   connected through a series of small number of nodes bonded covalently, i.e. L = {l1, l2, ... li, ... } 
   is a set of numInLoop nodes with li covalently bonded to both li-1 and li+1 ; 
   l1 is covalently bonded to r1 and lk is covalently bonded to t1.  
   Additionally their should be no hydrogen bonds among L (except possibly among neighbors). 
   Thus the set L forms the loop of the stem-loop {S1, S2, L}.

   Returns a sequence of stem-loops. Each stem-loop is a tuple. First element in a stem, second the loop.
  "
  [RAG len-stem num-in-loop]
  (filter (comp not nil?) 
    (map (fn [stem] (when-let [s-loop (stem-is-stem-loop RAG stem num-in-loop)]
                      [stem s-loop])) 
      (stems RAG len-stem))))

(def stem-loops (memoize stem-loops))

; stem-loops query
(comment
  (doseq [l (range 1 20)]
    (doseq [s (range 1 15)]
      (let [c (count (stem-loops RAG s l))]
        (if (> c 0) (println s l c))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;; KISSING LOOPS ;;;;;;;;;;;;;;;;;;;;;;;;

(defn kissing-loops
  "Returns all kissing-stem-loops.
   Definition of a kissing-stem-loop :
    * Two stem loops kiss when their loop sets interact via hydrogen bonds.
   Returns a 3-ple with first two elems the two stem loops, third elem the hydrogen bond forming nodes in the 2nd stem-loop"
  [RAG len-stem-1 num-in-loop-1 len-stem-2 num-in-loop-2]
  (let [stem-loops-1 (stem-loops RAG len-stem-1 num-in-loop-1)
        stem-loops-2 (stem-loops RAG len-stem-2 num-in-loop-2)]
    (filter (fn [k] (> (count (k 2)) 0)) 
      (for [sl-1 stem-loops-1, sl-2 stem-loops-2
            :when (let [t-sl-1 ((sl-1 0) 1)
                        s-sl-2 ((sl-2 0) 0)
                        max-t-sl-1 (t-sl-1 (dec (count t-sl-1)))
                        min-s-sl-2 (s-sl-2 (dec (count s-sl-2)))]
                    (< (res-id RAG max-t-sl-1) (res-id RAG min-s-sl-2)))]
        (let [sloop-sl-1 (sl-1 1)
              sloop-sl-2 (sl-2 1)
              hbonds-sl-1 (apply union (map #(hbond-edges RAG %) sloop-sl-1))]
          [sl-1 sl-2 (intersection hbonds-sl-1 (set sloop-sl-2))])))))

(def kissing-loops (memoize kissing-loops))

; kissing loops query
(comment
  (doseq [s1 (range 2 20)] 
    (doseq [l1 (range 1 20)] 
      (doseq [s2 (range s1 20)] 
        (doseq [l2 (range l1 20)] 
          (let [kl (kissing-loops RAG s1 l1 s2 l2)
                c (count kl)]
            (when (> c 0)
              (println s1 l1 s2 l2 c))))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;; BULGES ;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn two-stems-bulge
  "Checks if two stems have a buldge"
  [RAG s-1 s-2 bulge-len]
  (let [s-1-left (s-1 0)
        s-2-left (s-2 0)
        left-matches? (= ((cbond-edges RAG (s-2-left 0)) 1) (s-1-left (dec (count s-1-left))))]
    (when left-matches?
      (let [s-1-right (s-1 1)
            s-2-right (s-2 1)
            right-connection (connect RAG (s-1-right (dec (count s-1-right))) (s-2-right 0) bulge-len)]
        right-connection))))

(defn bulges
  "Returns all bulges.

   Definition of bulge:
   Bulges are two stems {{r1, r2, ...} {t1, t2, ...}} (length N1) and {{u1, u2, ...} {v1, v2, ...}} (length N2) 
    * such that rN1 is covalently bonded to u1 
    * and tN1 and v1 are on the same chain connected through a series of small number of nodes bonded covalently, 
      i.e. L = {l1, l2, ... li, ...} is a set of ‘K’ nodes with li covalently bonded to both li-1 and li+1 ; 
      l1 is covalently bonded to tN1 and lk is covalently bonded to v1. 
   Thus {{{r1, r2, ... } {t1, t2, ...}} {{u1, u2, ...} {v1, v2, ...}} {l1, l2, ... li, ...}} (N1, N2, K) is a bulge.

   Each bulge a 3-ple with first two elems the two stems, third elem is the bulge loop.
  "
  [RAG stem-len-1 stem-len-2 bulge-len]
  (let [stems-1 (stems RAG stem-len-1)
        stems-2 (stems RAG stem-len-2)]
    (filter (fn [b] ((comp not nil?) (b 2))) 
      (for [s-1 stems-1, s-2 stems-2]
        (vector s-1 s-2 (two-stems-bulge RAG s-1 s-2 bulge-len))))))

(def bulges (memoize bulges))

; buldges query
(comment
  (doseq [s1 (range 2 16)] 
    (doseq [s2 (range s1 16)] 
      (doseq [b (range 1 20)] 
        (let [bs (bulges RAG s1 s2 b)
              c (count bs)]
          (when (> c 0)
            (println s1 s2 b c)))))))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(comment
  (do
    (def RAG (form-RNA-graph {:atoms-pdb-file-name "C:\\Docs\\ramsagar\\1VQO-min.coor" 
                     :hbonds-file-name "C:\\Docs\\ramsagar\\1VQO-min.coor.out.grepped"}))
    (def RAG (form-RNA-graph {:atoms-pdb-file-name "C:\\Docs\\ramsagar\\Input.pdb" 
                                 :hbonds-file-name "C:\\Docs\\ramsagar\\Input.pdb.out.grepped"}))
    (get-RNA-graph-stats RAG)
    (print-RNA-graph-to-file RAG "test")
    (dotimes [i 34] (println (+ 2 i) "," (count (nlets-distinct RAG (+ 2 i)))))
    (dotimes [i 25] (println (+ 2 i) "," (count (stems RAG (+ 2 i)))))
    (doseq [l (range 1 20)]
      (doseq [s (range 1 15)]
        (let [c (count (stem-loops RAG s l))]
          (if (> c 0) (println s l c)))))
    (doseq [s1 (range 2 16)] 
      (doseq [s2 (range s1 16)] 
        (doseq [b (range 1 20)] 
          (let [bs (bulges RAG s1 s2 b)
                c (count bs)]
            (when (> c 0)
              (println s1 s2 b c))))))
    (doseq [s1 (range 2 20)] 
      (doseq [l1 (range 1 20)] 
        (doseq [s2 (range s1 20)] 
          (doseq [l2 (range l1 20)] 
            (let [kl (kissing-loops RAG s1 l1 s2 l2)
                  c (count kl)]
              (when (> c 0)
                (println s1 l1 s2 l2 c))))))))
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;; JAVA API Begins ;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn -formRNAGraph
  "Forms an instance of RNAGraph

   <pre>
   The RNAs graph is represented as (Adjacency list rep. - persistent) :
     1. nodes       : Map of Node names => Nodes
     2. cbond-edges : Map of Node names => Tuple of nodes [prev cbond, next cbond]
     3. hbond-edges : Map of Node names => Set of Node names

   Forms the RNA graph taking an options map with the following 2 keys -
   1. atoms-pdb-file-name
      -------------------
      Forms a sequence of residues from the Atoms PDB file. The information taken out is: 
      1.1 res-type (A, U, G, C) 
      1.2 res-num  (residue number in the pdb file)
      1.3 chain    (chain in the pdb file)
      Each residue has a unique name: R_&#60;chain&#62;&#60;res_num&#62;. Additionally each residue has a unique id.
      Residue nodes are 'Residue' objects having the following fields: name (R_&#60;chain&#62;&#60;res_num&#62;), type (A,U,G,C), chain, and id.
     
      After this the cbond edges are added to the graph. Residue with id k, is convalently bonded to residue with id k-1, k+1
      unless the residue is at the beginning/end of the chain.

   2. hbonds-file-name (Generated by Hbonds.java)
      ----------------
      Takes following info from the file - 
      2.1. Donor residue num, chain
      2.2. Acceptor residue num, chain
      *Note: do not change the format of the file as returned by Hbonds.java*
      Adds the hbonds edges to the graph as formed above.
   </pre>
   
   @param opts A hash map containing precisely the two keys - 'atoms-pdb-file-name', 'hbonds-file-name'

   @return     An immutable iiit.RNAGraph object

   @version    0.0.1
   @author     Sanchay Harneja"
  [opts]
  (form-RNA-graph opts))

(defn -getAllResidues
  "Get all residues from the RAG.

   @param RAG An instance of a RNAGraph

   @return    Returns a hash map of residue-names (R_&#60;chain&#62;&#60;res_num&#62;)=> Residue objects

   @version    0.0.1
   @author     Sanchay Harneja
  "
  [RAG]
  (:nodes RAG))

(defn -getAllCBonds
  "Get all CBond (undirected) edges from RAG. 

   @param RAG An instance of a RNAGraph

   @return    Returns a hash map of residue-names (R_&#60;chain&#62;&#60;res_num&#62;)=> Tuple of Residue names [prev-cbond, next-cbond]

   @version    0.0.1
   @author     Sanchay Harneja
  "
  [RAG]
  (:cbond-edges RAG))

(defn -getAllHBonds
  "Get all HBond (undirected) edges from RAG. 

   @param RAG An instance of a RNAGraph

   @return    Returns a hash map of residue-names (R_&#60;chain&#62;&#60;res_num&#62;) => Set of Residue names

   @version    0.0.1
   @author     Sanchay Harneja
  "
  [RAG]
  (:hbond-edges RAG))

(defn -getRNAGraphStats    
  "Get stats of an instance of RNAGraph
   
   <pre>
   Returns a map containing the following keys :
    - Number of Residue nodes
    - Number of Directed CBond edges
    - Number of Directed HBond edges
   </pre>

   @param RAG  An instance of a RNAGraph

   @return     A map containing the stats of the RNAGraph

   @version    0.0.1
   @author     Sanchay Harneja
   "
  [RAG]
  (get-RNA-graph-stats RAG))

(defn -printRNAGraphToFile
  "Print an instance of RNAGraph to a file
   
   <pre>
   Prints a RAG in a file named by fileName.  The file is overwritten every time.
   It prints in order :
    * Residue Nodes (1 per line)
    * CBond edges (1 residue's per line)
    * Hbond edges (1 residue's per line)
   </pre>

   @param RAG       An instance of a RNAGraph
   @param fileName  Name of the file to which to print the RAG. The file is overwritten every time.

   @version    0.0.1
   @author     Sanchay Harneja
   "
  [RAG fileName]
  (print-RNA-graph-to-file RAG fileName))

(defn -getNLets
  "Get all n-lets 

   <pre>
   Definition of nlet:
     A set of N connected (by hbonds) nodes.
   </pre>

   @param RAG  An instance of a RNAGraph
   @param n    The n in n-let

   @return     The Set of n-lets, each n-let is a sorted list of residues.

   @version    0.0.1
   @author     Sanchay Harneja
  "
  [RAG n]
  (nlets-distinct RAG n))

(defn -getStems
  "Get all stems of length stemLen.
   
   <pre>
   Definition of a stem of length n:
   A set S of two sets S1, S2 of equal number of 'n' nodes such that, 
    * S1 = {r1, r2, ... rn} with ri covalently bonded to both ri-1 and ri+1 ; 
    * S2 = {t1, t2, ... tn} with ti covalently bonded to both ti-1 and ti+1 ; 
    * Additionally each ri and ti for all i=1-to-n form hydrogen bonds.
   </pre>
  
   @param RAG     An instance of a RNAGraph
   @param stemLen All stems of length stemLen will be returned

   @return        Returns the sequence of stems. Each stem is a tuple. First elem is one side, second element the other side.

   @version    0.0.1
   @author     Sanchay Harneja
  " 
  [RAG stemLen]
  (stems RAG stemLen))

(defn -getStemLoops
  "Get stemLoops of length stemLen and loopLen residues in the loop. A.k.a hairPinLoops

   <pre>
   Definition of a stemLoop :
   Stem loops are stems of length stemLen with an additional properties 
    * r1 and t1 are on the same chain connected through a series of small number of nodes bonded covalently,
    * i.e. L = {l1, l2, ... li, ... } is a set of loopLen nodes with li covalently bonded to both li-1 and li+1 ; 
    * l1 is covalently bonded to r1 and lk is covalently bonded to t1.  
    * Additionally their should be no hydrogen bonds among L (except possibly among neighbors). 
   Thus the set L forms the loop of the stem-loop {S1, S2, L}.
   </pre>
   
   @param RAG     An instance of a RNAGraph
   @param stemLen Length of the stem
   @param loopLen Length of the loop

   @return        Returns the sequence of all stemsLoops. Each stem-loop is a tuple. First element in a stem, second the loop.

   @version    0.0.1
   @author     Sanchay Harneja
  "
  [RAG stemLen loopLen]
  (stem-loops RAG stemLen loopLen))


(defn -getBulges
  "Returns all bulges with first stem of length stem1Len, second of stem2Len and the bulge loop of len bulgeLen.

   <pre>
   Definition of bulge:
   Bulges are two stems {{r1, r2, ...} {t1, t2, ...}} (length stem1Len) and {{u1, u2, ...} {v1, v2, ...}} (length stem2Len) 
    * such that rN1 is covalently bonded to u1 
    * and tN1 and v1 are on the same chain connected through a series of small number of nodes bonded covalently, 
      i.e. L = {l1, l2, ... li, ...} is a set of bulgeLen nodes with li covalently bonded to both li-1 and li+1 ; 
      l1 is covalently bonded to tN1 and lk is covalently bonded to v1. 
   Thus {{{r1, r2, ... } {t1, t2, ...}} {{u1, u2, ...} {v1, v2, ...}} {l1, l2, ... li, ...}} (stem1Len, stem2Len, bulgeLen) is a bulge.
   </pre>

   @param RAG      An instance of a RNAGraph
   @param stem1Len Length of the first stem
   @param stem2Len Length of the second stem
   @param bulgeLen Length of the bulge loop

   @return        Returns the sequence of all bulges. Each bulge a 3-ple with first two elems the two stems, third elem is the bulge loop.

   @version    0.0.1
   @author     Sanchay Harneja

  "
  [RAG stem1Len stem2Len bulgeLen]
  (bulges RAG stem1Len stem2Len bulgeLen))

(defn -getKissingStemLoops 
  "Get all kissing stem loops.

   <pre>
   Definition of a kissing-stem-loop :
    * Two stem loops kiss when their loop sets interact via hydrogen bonds.
   </pre>

   @param RAG      An instance of a RNAGraph
   @param stem1Len Length of the stem in first  stem loop
   @param loop1Len Length of the loop in first  stem loop
   @param stem2Len Length of the stem in second stem loop
   @param loop2Len Length of the loop in second stem loop

   @return         Returns the sequence of all kissing stemsLoops. Each stem-loop is a 3-ple with first two elems the two stem loops, third elem the hydrogen bond forming nodes in the 2nd stem-loop.

   @version    0.0.1
   @author     Sanchay Harneja
   "
  [RAG stem1Len loop1Len stem2Len loop2Len]
  (kissing-loops RAG stem1Len loop1Len stem2Len loop2Len))


(defn- convert-fn-meta-into-javadoc
  [impl-ns prefix mname pclasses rclass is-static?]
  (let [m (meta (resolve (symbol (str (str impl-ns) "/" prefix (str mname)))))
        {:keys [file line doc]} m 
        arglist (first (filter #(= (count pclasses) (count %)) (:arglists m)))
        method-sign (str "public " 
                        (if is-static? "static " "") 
                        (str rclass) " "
                        (str mname) " ("
                        (apply str (drop-last (interleave pclasses (repeat " ") arglist (repeat ", "))))
                        ") {}")]
    (str
      "  /**\n"
      "   *  " (str-join "\n   *  " (re-split #"\n" doc)) "\n"
      "   *  Definition present at " {:file file :line line}  "\n"
      "   */\n  "
      method-sign
      "\n\n")))

(defn- generate-javadoc
  [options-map]
  (let [default-options {:prefix "-" :impl-ns (ns-name *ns*)}
        {:keys [name methods prefix impl-ns class-doc]}  (merge default-options options-map)
        [_ package c-name] (re-matches #"(.+)\.([^\.]+)" name)
        javadoc (str 
                  "package " package ";\n"
                  (str
                    "  /**\n"
                    "   *  " (str-join "\n   *  " (re-split #"\n" class-doc)) "\n"
                    "   */\n  "
                    "\n\n")
                  "public class " c-name "{\n\n")
        javadoc (loop [methods methods
                       javadoc javadoc]
                  (if-let [[mname pclasses rclass :as msig] (first methods)]
                    (recur (rest methods) 
                      (str javadoc
                        (convert-fn-meta-into-javadoc impl-ns prefix mname pclasses rclass (:static (meta msig)))))
                    javadoc))
        javadoc (str javadoc "}\n")]
    [package c-name javadoc]))

(defmacro gen-class+javadoc
  [& options]
  (do
    (macroexpand `(gen-class ~@options))
    (when *compile-files*
      (let [options-map              (into {} (map vec (partition 2 options)))
            [package c-name javadoc] (generate-javadoc options-map)
            path-of-java-file        (str (System/getProperty "user.dir") 
                                       java.io.File/separator "src" java.io.File/separator
                                       (.replaceAll package "\\." (str java.io.File/separator java.io.File/separator))
                                       java.io.File/separator c-name ".java")]
        (do
          (clojure.contrib.io/make-parents (java.io.File. path-of-java-file))
          (clojure.contrib.io/spit path-of-java-file javadoc)
          (let [cmd (str "javadoc -sourcepath src -d doc " package " > doc/javadoc.out 2> doc/javadoc.err")
                javadoc-ret (clojure.contrib.shell/sh "cmd.exe" "/c" cmd :return-map true)]
            ;(println javadoc-ret)
            )
          (clojure.contrib.io/delete-file path-of-java-file))))))

(println "Before gen-classing iiit.Topos")

(gen-class+javadoc
  :class-doc     "
                  This Class provides a purely stateless, persistent (in the immutable sense) API for - <br>
                  <ul>
                    <li> Forming a graph data structure for a RNA given its residues (nodes) and cbond, hbond (edges) information.</li>
                    <li> It provides means to extract information out the RAG (RNA as Graph) data structure thus formed.</li>
                    <li> Additionally it provides implementations to find certain Topological structures like n-lets, stems, stem-loops, bulges, kissing-stem-loops.</li>
                  </ul>
                  This API has been implemented in Clojure. <br><br>
                  For usage of this API see {@link iiit.ToposDemo iiit.ToposDemo}  <br><br>
                  For anything regarding this API, from understanding, documentation, examples, extension contact me sanchay.h@gmail.com
                 "
  :name          "iiit.Topos"
  :methods       [
                  #^{:static true} [formRNAGraph         [clojure.lang.PersistentArrayMap]  iiit.RNAAsGraph]
                  #^{:static true} [getAllResidues       [iiit.RNAAsGraph]                  clojure.lang.PersistentHashMap]
                  #^{:static true} [getAllCBonds         [iiit.RNAAsGraph]                  clojure.lang.PersistentHashMap]
                  #^{:static true} [getAllHBonds         [iiit.RNAAsGraph]                  clojure.lang.PersistentHashMap]
                  #^{:static true} [getRNAGraphStats     [iiit.RNAAsGraph]                  clojure.lang.PersistentArrayMap]
                  #^{:static true} [printRNAGraphToFile  [iiit.RNAAsGraph String]           void]
                  #^{:static true} [getNLets             [iiit.RNAAsGraph int]              clojure.lang.PersistentHashSet]
                  #^{:static true} [getStems             [iiit.RNAAsGraph int]              clojure.lang.LazySeq]
                  #^{:static true} [getStemLoops         [iiit.RNAAsGraph int int]          clojure.lang.LazySeq]
                  #^{:static true} [getBulges            [iiit.RNAAsGraph int int int]      clojure.lang.LazySeq]
                  #^{:static true} [getKissingStemLoops  [iiit.RNAAsGraph int int int int]  clojure.lang.LazySeq]
                 ]
  )

(println "Done loading iiit")






