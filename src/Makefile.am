AM_CPPFLAGS = -I$(top_builddir)/lib -I$(top_srcdir)/lib

bin_PROGRAMS = indel-seq-gen

indel_seq_gen_SOURCES = motif.h motif.cpp inTree.cpp inTree.h seqGenOptions.cpp seqGenOptions.h main.cpp evolve.cpp evolve.h model.cpp model.h nucmodels.cpp nucmodels.h aamodels.cpp aamodels.h eigen.cpp eigen.h treefile.cpp treefile.h tree.h progress.cpp progress.h gamma.cpp gamma.h twister.cpp twister.h random.cpp random.h insert_freqs.cpp insert_freqs.h trace.cpp trace.h inClade.cpp inClade.h prosite_motif.h
indel_seq_gen_LDADD = $(top_builddir)/lib/libgnu.a
