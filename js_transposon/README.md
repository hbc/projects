Performs merging and scoring from multiple mapped transposon experiments.

# Usage

Requires Java and commandline access to run. Download the latest release jar:
[hbc.transposon-0.0.2-standalone.jar][2] (25.3Mb).

Runs in a two step process. The first takes a YAML config or Excel sample file
and merges into a single output CSV file:

    $ java -jar hbc.transposon-0.0.2-standalone.jar 
      merge <work_directory> -c <YAML config file> -x <Excel sample file>
    
The second scores and filters the merged samples:
    
    $ java -jar hbc.transposon-0.0.2-standalone.jar 
      score <merged CSV file> -c <YAML config file> -x <Excel sample file>

# Development

Written in Clojure and the requires the [leiningen build tool][1] to run
directly from source:

    $ lein merge <work_directory> -c <YAML config file> -x <Excel sample file>
    $ lein score <merged CSV file> -c <YAML config file> -x <Excel sample file>

[1]: https://github.com/technomancy/leiningen
[2]: https://s3.amazonaws.com/hbc.transposon/hbc.transposon-0.0.2-standalone.jar
