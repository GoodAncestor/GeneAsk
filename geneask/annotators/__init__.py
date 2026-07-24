"""GeneAsk annotators — external variant-annotation sources, mirror-first.

Each annotator builds a local SQLite mirror from a public bulk download
(refresh()) and answers per-rsID lookups offline (get()). Same design as
MethylAsk's EWAS mirror. The heavy refresh runs on a worker node; the app
queries the built mirror.
"""
