CC     = gcc
CFLAGS = -std=c99 -g -I/usr/local/include/
EFLAGS =  
EFILE  = test_vbgmmfit
LIBS   = -lgomp -lpthread -lm -lgsl -lgslcblas -L/usr/local/lib
OBJS   = c_vbgmm_fit.o test_vbgmm_fit.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : c_vbgmm_fit.c c_vbgmm_fit.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o test_vbgmmfit
