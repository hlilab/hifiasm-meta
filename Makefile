CXX=		g++
CC=			gcc
CXXFLAGS=	-g3 -O3 -msse4.2 -mpopcnt -fomit-frame-pointer -Wall -Wno-unused -Wno-sign-compare
CFLAGS=		$(CXXFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		CommandLines.o Process_Read.o Assembly.o Hash_Table.o \
			POA.o Correct.o Levenshtein_distance.o Overlaps.o Trio.o kthread.o Purge_Dups.o \
			htab.o hist.o sketch.o anchor.o extract.o sys.o ksw2_extz2_sse.o \
			meta_util.o meta_util_debug.o Overlaps_hamt.o \
			sptree.o data.o tsne.o t-sne.o main.o
EXE=		hifiasm_meta
LIBS=		-lz -lpthread -lm

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.cpp .c .o
.PHONY:all clean depend

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

$(EXE):$(OBJS) main.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(EXE) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CPPFLAGS) $(DFLAGS) -- *.cpp)

# DO NOT DELETE

Assembly.o: Assembly.h CommandLines.h Process_Read.h Overlaps.h kvec.h kdq.h
Assembly.o: Hash_Table.h htab.h meta_util.h POA.h Correct.h
Assembly.o: Levenshtein_distance.h kthread.h ksort.h meta_util_debug.h
CommandLines.o: CommandLines.h ketopt.h gitcommit.h
Correct.o: Correct.h Hash_Table.h htab.h Process_Read.h Overlaps.h kvec.h
Correct.o: kdq.h CommandLines.h meta_util.h Levenshtein_distance.h POA.h
Correct.o: Assembly.h ksw2.h ksort.h
Hash_Table.o: Hash_Table.h htab.h Process_Read.h Overlaps.h kvec.h kdq.h
Hash_Table.o: CommandLines.h meta_util.h ksort.h
Levenshtein_distance.o: Levenshtein_distance.h
Output.o: Output.h CommandLines.h
Overlaps.o: Overlaps.h kvec.h kdq.h ksort.h Process_Read.h CommandLines.h
Overlaps.o: Hash_Table.h htab.h meta_util.h Correct.h Levenshtein_distance.h
Overlaps.o: POA.h Purge_Dups.h Overlaps_hamt.h kthread.h
Overlaps_hamt.o: Overlaps.h kvec.h kdq.h Process_Read.h CommandLines.h
Overlaps_hamt.o: Overlaps_hamt.h ksort.h htab.h kvec.h t-sne.h
POA.o: POA.h Hash_Table.h htab.h Process_Read.h Overlaps.h kvec.h kdq.h
POA.o: CommandLines.h meta_util.h Correct.h Levenshtein_distance.h
Process_Read.o: Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
Purge_Dups.o: ksort.h Purge_Dups.h kvec.h kdq.h Overlaps.h Hash_Table.h
Purge_Dups.o: htab.h Process_Read.h CommandLines.h meta_util.h Correct.h
Purge_Dups.o: Levenshtein_distance.h POA.h kthread.h
Trio.o: khashl.h kthread.h kseq.h Process_Read.h Overlaps.h kvec.h kdq.h
Trio.o: CommandLines.h htab.h meta_util.h
anchor.o: htab.h Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
anchor.o: meta_util.h ksort.h Hash_Table.h
extract.o: Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h khashl.h
extract.o: kseq.h
hist.o: htab.h Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
hist.o: meta_util.h
htab.o: kthread.h khashl.h kseq.h ksort.h htab.h Process_Read.h Overlaps.h
htab.o: kvec.h kdq.h CommandLines.h meta_util.h
kthread.o: kthread.h
main.o: CommandLines.h Process_Read.h Overlaps.h kvec.h kdq.h Assembly.h
main.o: Levenshtein_distance.h htab.h meta_util.h meta_util_debug.h
meta_util.o: meta_util.h Process_Read.h Overlaps.h kvec.h kdq.h
meta_util.o: CommandLines.h kthread.h
meta_util_debug.o: Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
sketch.o: kvec.h htab.h Process_Read.h Overlaps.h kdq.h CommandLines.h
sketch.o: meta_util.h
sys.o: htab.h Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
sys.o: meta_util.h
# bhtsne
data.o: t-sne.h kseq.h
