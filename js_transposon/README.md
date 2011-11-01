Performs merging and scoring from multiple mapped transposon experiments.

Written in Clojure and requires Java and the [lein build tool][1].

To run:

    $ lein deps
    $ lein run :merge <work_directory> <YAML config file>
    $ lein run :score <merged CSV file>

[1]: https://github.com/technomancy/leiningen
