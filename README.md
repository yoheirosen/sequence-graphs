##Dependencies

###To install Apache Spark with GraphX

```
git clone https://github.com/apache/incubator-spark
cd incubator-spark
git checkout v0.9.0-incubating
sbt publish-local
```

###To install vcfimp

```
git clone https://github.com/adamnovak/vcfimp.git
cd vcfimp
git checkout v0.6.1
sbt publish-local
```

##Installation

###Building

```
sbt stage
```

###Testing

```
sbt test
```

###Running command-line tools

```
./importVCF.sh <vcf file> <sample name>
```
