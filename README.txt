Modelled PDB data in Clojure to search for higher level topological patterns in the secondary structure of RNA.

COMPILE WINDOWS
---------------

..\IIIT_RNA_Project>java -cp lib\clojure-1.2.0-master-20100518.110252-71.jar;lib\clojure-contrib-1.2.0-20100518.110549-109.jar;lib\commons-lang-2.3.jar;classes;src clojure.main -e "(compile 'iiit)"

..\IIIT_RNA_Project>javac -classpath lib\clojure-1.2.0-master-20100518.110252-71.jar;lib\clojure-contrib-1.2.0-20100518.110549-109.jar;lib\commons-lang-2.3.jar;classes -d classes src\iiit\ToposDemo.java

..\IIIT_RNA_Project>javac -classpath lib\clojure-1.2.0-master-20100518.110252-71.jar;lib\clojure-contrib-1.2.0-20100518.110549-109.jar;lib\commons-lang-2.3.jar;classes -d classes src\iiit\HBonds.java

RUN WINDOWS
-----------

..\IIIT_RNA_Project>java -classpath lib\clojure-1.2.0-master-20100518.110252-71.jar;lib\clojure-contrib-1.2.0-20100518.110549-109.jar;lib\commons-lang-2.3.jar;classes iiit.ToposDemo

..\IIIT_RNA_Project>java -classpath lib\clojure-1.2.0-master-20100518.110252-71.jar;lib\clojure-contrib-1.2.0-20100518.110549-109.jar;lib\commons-lang-2.3.jar;classes iiit.HBonds


COMPILE LINUX
--------------

../IIIT_RNA_Project>java -cp lib/clojure-1.2.0-master-20100518.110252-71.jar:lib/clojure-contrib-1.2.0-20100518.110549-109.jar:lib/commons-lang-2.3.jar:classes:src clojure.main -e "(compile 'iiit)"

../IIIT_RNA_Project>javac -classpath lib/clojure-1.2.0-master-20100518.110252-71.jar:lib/clojure-contrib-1.2.0-20100518.110549-109.jar:lib/commons-lang-2.3.jar:classes -d classes src/iiit/ToposDemo.java

../IIIT_RNA_Project>javac -classpath lib/clojure-1.2.0-master-20100518.110252-71.jar:lib/clojure-contrib-1.2.0-20100518.110549-109.jar:lib/commons-lang-2.3.jar:classes -d classes src/iiit/HBonds.java

RUN LINUX
-----------

../IIIT_RNA_Project>java -classpath lib/clojure-1.2.0-master-20100518.110252-71.jar:lib/clojure-contrib-1.2.0-20100518.110549-109.jar:lib/commons-lang-2.3.jar:classes iiit.ToposDemo

../IIIT_RNA_Project>java -classpath lib/clojure-1.2.0-master-20100518.110252-71.jar:lib/clojure-contrib-1.2.0-20100518.110549-109.jar:lib/commons-lang-2.3.jar:classes iiit.HBonds
