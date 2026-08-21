#include "OEuDE.h"

//double oeude(double *solution, uint64_t cr, uint64_t budget, uint64_t *seed, int dim, double (*feval)(double *, int)){
double oeude(double *elite, uint64_t budget, int dim, double (*feval)(double *, int), uint64_t CR){
    uint64_t k, l, m, elite=0, best;
    int32_t i, j;
    double *pop, *fpop, F, bestValue;
    uint64_t ind=0x543210, unused_elements=0xba9876, copy_ind;
    uint64_t i_smp[6];
    uint64_t n_feval    = 0;
    fpop    = (double*)calloc(2*NPOP_OE, sizeof(double));
    if(fpop==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    pop     = (double*)calloc((2*NPOP_OE*dim), sizeof(double));
    if(pop==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    //init pop
    for(i=0; i<NPOP_OE; i++){
        for(j=0; j<dim; j++){
            pop[i*dim + j] = 2.0f*doubleRand() - 1.0f;
        }
        fpop[i]=feval(&pop[i*dim], dim); 
        n_feval++;
        if(i>0){//Sort population according to fitness values
            k   = i;
            for(j=i-1; j>=0; j--){
                if(fpop[i] < fpop[get_ind(ind, j)]){
                    k=j;
                }else{
                    break;
                }
            }
            if(i != k){
                ind = insert_ind(ind, k, i);
            }else{
                ind = insert_ind(ind, i, i);
            }
        }
    }

    do{
        //Oposition population calculation
        for(i=NPOP_OE-1; i>=0; i--){
            k   = get_ind(unused_elements, i);
            l   = get_ind(ind, i);
            for(j=0; j<dim; j++){
                pop[k*dim + j] = -pop[l*dim + j];
            }
            fpop[k] = feval(&pop[k*dim], dim); 
            n_feval++;
            if(budget <= n_feval){
                break;
            }
        }

        //Select the best elements
        for(i=NPOP_OE-1; i>=0; i--){
            k   = get_ind(unused_elements, i);
            m   = 9999; //if changed an update will be performed
            for(j=NPOP_OE-1; j>=0; j--){
                l   = get_ind(ind, j);
                if(fpop[k] < fpop[l]){
                    m=j;
                }else{
                    break;
                }
            }
            if(m != 9999){ //Check if there are modifications and sort population acordingly for the next generation
                ind             = insert_ind(ind, m, k);
                unused_elements = replace_ind(unused_elements, i, get_ind(ind, NPOP_OE));
            }
        }
        
        //DE over each population element
        best = get_ind(ind, 0)*dim;
        copy_ind    = ind;
        for(i=0; i < NPOP_OE; i++){
            k   = get_ind(unused_elements, 0);
            i_smp[0]    = i;
            randSample(i_smp, NPOP_OE, 5);
            j   = (next() % (uint64_t)dim);
            uint64_t trial = get_ind(unused_elements, 0)*dim;
            for (l=0; l<6; l++){
                i_smp[l] = get_ind(ind, i_smp[l])*dim;
            }
            for(l=0; l<dim; l++){
                if(next() < CR || l==j){
                    F   = 1.4f * doubleRand() + 0.1f;
                    m   = next() % 5;
                    switch (m) //Select operator on each component randomly
                    {
                    case 0://rand/1
                        pop[trial + l]  = pop[i_smp[1]+l] + F * (pop[i_smp[2]+l] - pop[i_smp[3]+l]);
                        break;
                    case 1://best/1
                        pop[trial + l]  = pop[best+l] + F * (pop[i_smp[1]+l] - pop[i_smp[2]+l]);
                        break;
                    case 2://target-to-best/1
                        pop[trial + l]  = pop[i_smp[0]+l] + F*(pop[best+l] - pop[i_smp[0]+l] + pop[i_smp[1]+l] - pop[i_smp[2]+l]);
                        break;
                    case 3://rand/2
                        pop[trial + l]  = pop[i_smp[1]+l] + F * (pop[i_smp[2]+l] - pop[i_smp[3]+l] + pop[i_smp[4]+l] - pop[i_smp[5]+l] );
                        break;
                    default://best/2
                        pop[trial + l]  = pop[best+l] + F * (pop[i_smp[1]+l] - pop[i_smp[2]+l] + pop[i_smp[3]+l] - pop[i_smp[4]+l] );
                        break;
                    }
                    //repair boundaries
                    while (pop[trial + l] > 1.0f)
                    {
                        pop[trial + l]-=2.0f;
                    }
                    while (pop[trial + l] < -1.0f)
                    {
                        pop[trial + l]+=2.0f;
                    }
                    
                }else{
                    pop[trial+l] = pop[i_smp[0]+l];
                }
            }
            fpop[k]=feval(&pop[trial], dim);
            n_feval++;
            //in case an element was outperformed is replaced and is accomodated in its final position through insert sorting
            if(fpop[k] < fpop[get_ind(ind, i)]){
                unused_elements = replace_ind(unused_elements, 0, get_ind(ind, i));
                ind             = replace_ind(ind, i, k);
                copy_ind        = replace_ind(copy_ind, i, k);
                if(i>0){
                    k   = i;
                    for(j=i-1; j>=0; j--){
                        if(fpop[get_ind(copy_ind, i)] < fpop[get_ind(copy_ind, j)]){
                            k=j;
                        }else{
                            break;
                        }
                    }
                    if(i != k){
                        copy_ind = insert_ind(copy_ind, k, get_ind(copy_ind, i));
                    }
                }
            }
            
            if(budget <= n_feval){
                break;
            }
        }
        ind = copy_ind;
    } while(budget > n_feval);
    //Copy best solution to the final array
    bestValue   = fpop[get_ind(ind, 0)];
    for(i=0; i<dim; i++){
        elite[i] = pop[get_ind(ind,0)*dim + i];
    }
    cleanup:
    free(pop);
    free(fpop);
    return bestValue;
}

int compareDouble(const void *a, const void *b) {
    double x = *(const double *)a;
    double y = *(const double *)b;

    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

uint64_t insert_ind(uint64_t ind, uint64_t pos, uint64_t value){
    uint64_t apos, find;
    apos   = (7 - pos) << 2;
    find   = (MAXINTB >> (apos)) & (ind);
    apos   = pos << 2;
    find   += (MAXINTA << (apos)) & (ind << 4);
    find   += (value<<apos);
    return find;
}

uint64_t replace_ind(uint64_t ind, uint64_t pos, uint64_t value){
    uint64_t apos, find;
    apos   = (7 - pos) << 2;
    find   = (MAXINTB >> (apos)) & (ind);
    apos   = pos << 2;
    find   += (MAXINTA << (apos)) & (ind);
    find   += (value<<apos);
    return find;
}

uint64_t move_ind(uint64_t ind, uint64_t pos0, uint64_t pos1){
    uint64_t apos, find, mask;
    apos    = pos0 << 2;
    mask    = MAXINTA << apos;
    apos    = (7 - pos1) << 2;
    mask    += (MAXINTB >> (apos));
    find    = ind & mask;
    mask    = ~mask;
    find    += ((mask & ind)<<4) & mask;
    find    += (((mask & ind) >> ((pos0 - pos1)<<2))& mask);
    return find;

}

uint64_t get_ind(uint64_t ind, uint64_t pos){
    uint64_t ret;
    ret = (pos<<2);
    ret = (ind >> ret) & MASKPOS;
    return ret;
}