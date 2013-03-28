Help on method translate in module Bio.Seq:

translate(self, table='Standard', stop_symbol='*', to_stop=False, cds=False) method of Bio.Seq.Seq instance
    Turns a nucleotide sequence into a protein sequence. New Seq object.

    This method will translate DNA or RNA sequences, and those with a
    nucleotide or generic alphabet.  Trying to translate a protein
    sequence raises an exception.

    Arguments:
     - table - Which codon table to use?  This can be either a name
               (string), an NCBI identifier (integer), or a CodonTable
               object (useful for non-standard genetic codes).  This
               defaults to the "Standard" table.
     - stop_symbol - Single character string, what to use for terminators.
                     This defaults to the asterisk, "*".
     - to_stop - Boolean, defaults to False meaning do a full translation
                 continuing on past any stop codons (translated as the
                 specified stop_symbol).  If True, translation is
                 terminated at the first in frame stop codon (and the
                 stop_symbol is not appended to the returned protein
                 sequence).
     - cds - Boolean, indicates this is a complete CDS.  If True,
             this checks the sequence starts with a valid alternative start
             codon (which will be translated as methionine, M), that the
             sequence length is a multiple of three, and that there is a
             single in frame stop codon at the end (this will be excluded
             from the protein sequence, regardless of the to_stop option).
             If these tests fail, an exception is raised.

    e.g. Using the standard table:

    >>> coding_dna = Seq("GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    >>> coding_dna.translate()
    Seq('VAIVMGR*KGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
    >>> coding_dna.translate(stop_symbol="@")
    Seq('VAIVMGR@KGAR@', HasStopCodon(ExtendedIUPACProtein(), '@'))
    >>> coding_dna.translate(to_stop=True)
    Seq('VAIVMGR', ExtendedIUPACProtein())

    Now using NCBI table 2, where TGA is not a stop codon:

    >>> coding_dna.translate(table=2)
    Seq('VAIVMGRWKGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
    >>> coding_dna.translate(stop_symbol="@")
    Seq('VAIVMGR@KGAR@', HasStopCodon(ExtendedIUPACProtein(), '@'))
    >>> coding_dna.translate(to_stop=True)
    Seq('VAIVMGR', ExtendedIUPACProtein())

    Now using NCBI table 2, where TGA is not a stop codon:

    >>> coding_dna.translate(table=2)
    Seq('VAIVMGRWKGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
    >>> coding_dna.translate(table=2, to_stop=True)
    Seq('VAIVMGRWKGAR', ExtendedIUPACProtein())

    In fact, GTG is an alternative start codon under NCBI table 2, meaning
    this sequence could be a complete CDS:

    >>> coding_dna.translate(table=2, cds=True)
    Seq('MAIVMGRWKGAR', ExtendedIUPACProtein())

    It isn't a valid CDS under NCBI table 1, due to both the start codon and
    also the in frame stop codons:

    >>> coding_dna.translate(table=1, cds=True)
    Traceback (most recent call last):
        ...
    TranslationError: First codon 'GTG' is not a start codon

    If the sequence has no in-frame stop codon, then the to_stop argument
    has no effect:

    >>> coding_dna2 = Seq("TTGGCCATTGTAATGGGCCGC")
    >>> coding_dna2.translate()
    Seq('LAIVMGR', ExtendedIUPACProtein())
    >>> coding_dna2.translate(to_stop=True)
    Seq('LAIVMGR', ExtendedIUPACProtein())

    NOTE - Ambiguous codons like "TAN" or "NNN" could be an amino acid
    or a stop codon.  These are translated as "X".  Any invalid codon
    (e.g. "TA?" or "T-A") will throw a TranslationError.

    NOTE - Does NOT support gapped sequences.

    NOTE - This does NOT behave like the python string's translate
    method.  For that use str(my_seq).translate(...) instead.

