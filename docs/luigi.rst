
Comparison with luigi
=====================

The major source of inspiration for the `apetype` package was `luigi`
developed at `Spotify`. When it was developed, type hints were not a
thing yet in `Python`, so it is very understandable that it was not
used in `luigi`. Using a lot of luigi.Parameters makes the code look
bloated and the way of inheritance between tasks did also not give
issues. With that in mind, `apetype` was developed. In this section, a
1 on 1 comparison with the toy example of `luigi`.

Top artists example
-------------------

The code snippets on luigi's documentation page
( https://luigi.readthedocs.io/en/stable/example_top_artists.html )
are not directly executable, but the full code can be found at
``examples/top_artists.py`` in the `luigi` repository. In the
comparison here, code should also not be executed. It is just to
compare the syntax.


Aggregate Artist Streams with `luigi`
*************************************

    >>> class AggregateArtists(luigi.Task):
    ...     date_interval = luigi.DateIntervalParameter()
    ... 
    ...     def output(self):
    ...         return luigi.LocalTarget("data/artist_streams_%s.tsv" % self.date_interval)
    ... 
    ...     def requires(self):
    ...         return [Streams(date) for date in self.date_interval]
    ... 
    ...     def run(self):
    ...         artist_count = defaultdict(int)
    ... 
    ...         for input in self.input():
    ...             with input.open('r') as in_file:
    ...                 for line in in_file:
    ...                     timestamp, artist, track = line.strip().split()
    ...                     artist_count[artist] += 1
    ... 
    ...         with self.output().open('w') as out_file:
    ...             for artist, count in artist_count.iteritems():
    ...                 print(artist, count, file=out_file)


Aggregate Artist Streams with `apetype`
***************************************

    >>> import apetype as at
    ... import datetime, pathlib
    ... 
    ... class DateInterval(at.TaskBase):
    ...     days: float = 0
    ...     hours: float = 0
    ... 
    ...     def timedelta(_) -> datetime.timedelta:
    ...         return datetime.timedelta(days=_.days, hours=_.hours)
    ...
    ... class Streams(at.TaskBase):
    ...     date: float
    ... 
    ... class AggregateArtists(at.TaskBase):
    ...     date_interval: DateInterval
    ...
    ...     def dates_list(_, date_interval) -> list:
    ...         return [d for d in date_interval.timedelta]
    ...  
    ...     def output(_, date_interval) -> pathlib.Path:
    ...         return pathlib.Path("data/artist_streams_%s.tsv" % date_interval)
    ... 
    ...     def streams_list(_, dates_list: at.tasks.InjectItems) -> list:
    ...         return Streams(dates_list)
    ... 
    ...     def main(_, streams_list, output):
    ...         artist_count = defaultdict(int)
    ... 
    ...         for input in streams_list:
    ...             with input.open('r') as in_file:
    ...                 for line in in_file:
    ...                     timestamp, artist, track = line.strip().split()
    ...                     artist_count[artist] += 1
    ... 
    ...         with output.open('w') as out_file:
    ...             for artist, count in artist_count.iteritems():
    ...                 print(artist, count, file=out_file)



Conclusion
----------

Watch out Luigi, Donkey Kong is back and he wants to take over the
plumbing business .. in Python programming, with a keen sense and
(in)appropriate use of class.
