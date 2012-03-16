AM_CPPFLAGS = -I$(top_builddir)/lib -I$(top_srcdir)/lib

bin_PROGRAMS = indel-seq-gen

indel_seq_gen_SOURCES = motif.h motif.cpp inTree.cpp inTree.h seqGenOptions.cpp seqGenOptions.h main.cpp evolve.cpp evolve.h model.cpp model.h eigen.cpp eigen.h treefile.cpp treefile.h tree.h gamma.cpp gamma.h random.cpp random.h insert_freqs.cpp insert_freqs.h trace.cpp trace.h inClade.cpp inClade.h prosite_motif.h paleo.h paleo.cpp propose_path.cpp propose_path.h stats.cpp stats.h dependency.h dependency.cpp forward_simulation.cpp forward_simulation.h twister.h twister.cpp
indel_seq_gen_LDADD = $(top_builddir)/lib/libgnu.a
