Performs merging and scoring from multiple mapped transposon experiments.

# Usage

Written in Clojure and requires Java and the [lein build tool][1].

Runs in a two step process. The first takes a YAML config or Excel sample file
and merges into a single output CSV file:

    $ lein merge <work_directory> -c <YAML config file> -x <Excel sample file>
    
The second scores and filters the merged samples:
    
    $ lein score <merged CSV file> -c <YAML config file> -x <Excel sample file>

[1]: https://github.com/technomancy/leiningen
