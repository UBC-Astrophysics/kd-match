CFLAGS = -g
GCC = gcc
.c.o :
	$(GCC) -c $(CFLAGS) $*.c
EXES = match_kd pair_kd triangle_kd quad_kd calctrans transform
all : $(EXES)
MATCHOBJS =  match_kd.o kdtree.o
match_kd : $(MATCHOBJS) 
	gcc $(CFLAGS) -o match_kd $(MATCHOBJS) -lm 
PAIROBJS = pair_kd.o kdtree.o loadfile.o
pair_kd : $(PAIROBJS)
	gcc $(CFLAGS) -o pair_kd $(PAIROBJS) -lm 
TRIOBJS = triangle_kd.o calctransform.o kdtree.o loadfile.o
triangle_kd : $(TRIOBJS)
	gcc $(CFLAGS) -o triangle_kd $(TRIOBJS) -lm 
QUADOBJS = quad_kd.o calctransform.o kdtree.o loadfile.o
quad_kd : $(QUADOBJS)  
	gcc $(CFLAGS) -o quad_kd $(QUADOBJS) -lm 
CALCOBJS = calctrans.o loadfile.o calctransform.o
calctrans : $(CALCOBJS) 
	gcc $(CFLAGS) -o calctrans $(CALCOBJS) -lm 
TRANSFORMOBJS = transform.o 
transform : $(TRANSFORMOBJS)
	gcc $(CFLAGS) -o transform $(TRANSFORMOBJS) -lm 
clean :
	rm *.o
dist-clean : clean
	rm $(EXES)
package.tar.gz : *.[ch] README Makefile
	tar cf package.tar *.[ch] README Makefile
	gzip -f package.tar
