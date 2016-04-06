/* used by heap.c */

#ifndef NHEAP_HH
#define NHEAP_HH

#define MAXELEMS 1000

#include <iostream>
#include <vector>

using namespace std;

template <class V>
class fHeap {
public:
    typedef struct el
    {
        V value;
        long idx;
    } elem;

    struct elemcomp {
        bool operator () (elem left, elem right) const
        { return left.value<right.value; }
    };

    enum sign {INC=0, DEC=1};

    /*returns the parent of a heap position*/
//    inline int parent(int posel) {
//        if (posel%2) return posel/2;
//        else return (posel-1)/2;
//    };

    inline long parent(long posel) {
        return (posel-1) >> 1;
    }

    fHeap(int sign=INC) {
        num_elems = 0;
        this->sign=sign;
    };

    ~fHeap() {
        vector<elem>().swap(heap);
        vector<long>().swap(order);
    };

    void clear() {
        heap.clear();
        order.clear();
        num_elems = 0;
    }

    bool isExisted(long idx, V &v) {
        if (num_elems>0 && idx<order.size() && order[idx]!=-1) {
            v=heap[order[idx]].value;
            return true;
        } else
            return false;
    };

    bool isExisted(long idx) {
        if (num_elems>0 && idx<order.size() && order[idx]!=-1) {
            return true;
        } else
            return false;
    };

    // old remove function with bug!!!

    //bool remove(int idx) {
    //	if (num_elems>0 && idx<order.size() && order[idx]!=-1) {
    //		int posel = order[idx];
    //           heap[posel].value = -1;
    //           moveup(posel);

    //           V value;
    //           dequeue(idx, value);

    //		return true;
    //	}
    //	return false;
    //};


    bool remove(long idx)
    {
        if (num_elems > 0 && idx < order.size() && order[idx] != -1)
        {
            long pos = order[idx];
            order[idx] = -1;
            if (num_elems - 1 > 0)
            {
                order[heap[num_elems-1].idx] = pos;
                heap[pos].idx = heap[num_elems-1].idx;
                //heap[pos].value = heap[endpos].value;
                num_elems--;
                updatequeue(heap[pos].idx, heap[num_elems].value);
            }
            else num_elems--;
            return true;
        }
        return false;
    }

    void enqueue(long idx, V value) {
        //if (num_elems>0 && idx<order.size() && order[idx]!=-1) {
        //    updatequeue(idx, value);
        //    return;
        //}

        V tmp;
        long tmpidx,p,posel;

        posel = num_elems; //last position
        if (heap.size()<num_elems+1)
            heap.push_back(elem());
        heap[num_elems].value = value;
        heap[num_elems++].idx = idx;
        if (idx>=order.size())
            order.insert(order.end(), idx-order.size()+1, -1);
        order[idx] = posel;

        while (posel >0) {
            p=parent(posel);
            //if (value<heap[p].value) {
            if (compare(value, heap[p].value)) {
                /* swap element with its parent */
                tmp = heap[p].value;
                heap[p].value = value;
                heap[posel].value = tmp;

                tmpidx = heap[p].idx;
                heap[p].idx = idx;
                heap[posel].idx = tmpidx;

                /* adjust order */
                order[idx] = p;
                order[tmpidx] = posel;

                posel = parent(posel);
            } else
                break;
        }
    };

    long dequeue(long &idx) {
        if (num_elems==0) /* empty queue */
            return 0;

        //value = heap[0].value;
        idx = heap[0].idx;
        //printf("dequeue %d", idx);
        order[idx] = -1; /* not in queue */
        if (num_elems-1>0) {
            heap[0].value = heap[num_elems-1].value;
            heap[0].idx = heap[num_elems-1].idx;
            order[heap[0].idx] = 0;
            num_elems--;
            movedown(0);
        } else
            num_elems--;

        return 1;
    };

    long dequeue() {
        if (num_elems==0) /* empty queue */
            return 0;

        //value = heap[0].value;
        long idx = heap[0].idx;
        //printf("dequeue %d ", idx);
        order[idx] = -1; /* not in queue */
        if (num_elems-1>0) {
            heap[0].value = heap[num_elems-1].value;
            heap[0].idx = heap[num_elems-1].idx;
            order[heap[0].idx] = 0;
            num_elems--;
            movedown(0);
        } else
            num_elems--;

        return 1;
    };

    int dequeue(long &idx, V &value) {
        if (num_elems==0) /* empty queue */
            return 0;

        value = heap[0].value;
        idx = heap[0].idx;
        order[idx] = -1; /* not in queue */
        if (num_elems-1>0) {
            heap[0].value = heap[num_elems-1].value;
            heap[0].idx = heap[num_elems-1].idx;
            order[heap[0].idx] = 0;
            num_elems--;
            movedown(0);
        } else
            num_elems--;

        return 1;
    };

    void updatequeue(long idx, V new_val) {
        long posel;

        posel = order[idx];

        if (new_val==heap[posel].value)
            return;

        //if (new_val<=heap[posel].value)
        if (!compare(heap[posel].value, new_val))
        {
            heap[posel].value = new_val;
            moveup(posel);
        } else
        {
            heap[posel].value = new_val;
            movedown(posel);
        }
    };

    void movedown(long md_start) {
        V tmp;
        long tmpidx;
        long posel = md_start; //root
        long swap;
        /*while posel is not a leaf and heap[posel].value > any of childen*/
        while ((posel<<1)+1 < num_elems) { // exists a left son
            if ((posel<<1)+2 < num_elems) { // exists a right son
                //if (heap[posel*2+1].value<heap[posel*2+2].value) // find minimum son
                if (compare(heap[(posel<<1)+1].value, heap[(posel<<1)+2].value))
                    swap = (posel<<1)+1;
                else
                    swap = (posel<<1)+2;
            } else
                swap = (posel<<1)+1;

            //if (heap[posel].value > heap[swap].value) { /*larger than smallest son*/
            if (compare(heap[swap].value, heap[posel].value)) {
                /*swap elements*/
                tmp = heap[swap].value;
                heap[swap].value = heap[posel].value;
                heap[posel].value = tmp;

                tmpidx = heap[swap].idx;
                heap[swap].idx = heap[posel].idx;
                heap[posel].idx = tmpidx;

                /*swap orders*/
                order[heap[swap].idx] = swap;
                order[heap[posel].idx] = posel;

                posel = swap;
            } else {
                break;
            }
        }
    };

    void moveup(long md_start=0) {
        V tmp;
        long tmpidx;
        long posel = md_start; //root
        long p;
        long idx=heap[posel].idx;

        /*while posel is not a root and heap[posel].value < parent */
        while (posel >0) {
            p=parent(posel);
            //if (heap[posel].value<heap[p].value) {
            if (compare(heap[posel].value,heap[p].value)) {
                /* swap element with its parent */
                tmp = heap[p].value;
                heap[p].value = heap[posel].value;
                heap[posel].value = tmp;

                tmpidx = heap[p].idx;
                heap[p].idx = heap[posel].idx;
                heap[posel].idx = tmpidx;

                /* adjust order */
                order[heap[p].idx] = p;
                order[heap[posel].idx] = posel;

                posel = p;
            } else
                break;
        }
    };

    void getTop(long &idx, V &value) {
        idx=heap[0].idx;
        value=heap[0].value;
    };

    void reset() {
        heap.clear();
        order.clear();
        num_elems=0;
    };

    void reset2() {
        heap.clear();
        num_elems=0;
    };

    V getTopValue() { return heap[0].value; };
    long getTopIdx() { return heap[0].idx; };
    V getSecondTopValue() {
        return heap[1].value < heap[2].value ? heap[1].value : heap[2].value;
    }

    long size() { return num_elems; }

    void prlong_heap() {
        for (long i=0; i<num_elems; i++)
            printf("(%d) ", order[heap[i].idx]);
        printf("\n");
        //system("pause");
    };

private:
    bool compare(V a, V b) { return (sign==INC)? (a<b) : (a>b); };

    vector<elem> heap;
    vector<long> order;
    long num_elems;
    long sign;
};

///
/// Heap - end

#endif
