                                       SnpEff
                                       ------
                                       
Documentation

	http://snpeff.sourceforge.net/

Enjoy!

	Pablo Cingolani


To build SnpEff for Spritz:
1. Activate spritz conda environment
2. `conda install maven`
3. Execute the following commands:

```
cd SnpEff
mvn install:install-file -Dfile=lib/antlr-4.5.1-complete.jar -DgroupId=org.antlr -DartifactId=antlr -Dversion=4.5.1 -Dpackaging=jar
mvn install:install-file -Dfile=lib/biojava3-core-3.0.7.jar -DgroupId=org.biojava -DartifactId=biojava3-core -Dversion=3.0.7 -Dpackaging=jar
mvn install:install-file -Dfile=lib/biojava3-structure-3.0.7.jar -DgroupId=org.biojava -DartifactId=biojava3-structure -Dversion=3.0.7 -Dpackaging=jar
export VERSION=4.3
export VERSION_UND=`echo $VERSION | tr '.' '_'`
mvn clean compile assembly:assembly
mvn install:install-file -Dfile=target/SnpEff-$VERSION.jar -DgroupId=org.snpeff -DartifactId=SnpEff -Dversion=$VERSION -Dpackaging=jar -DgeneratePom=true --quiet
cp target/SnpEff-$VERSION-jar-with-dependencies.jar snpEff.jar
cd ..
```

4. To create new zip-file release, include everything except the `data` folder.
