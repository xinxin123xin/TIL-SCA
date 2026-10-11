/* Replay archived observations through the paper's extracted accumulators. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/reducers.h"

typedef struct {
    char name[96];
    Stats hw[2][2];
    Mixed pair[2];
} Campaign;

int main(int argc, char **argv) {
    if (argc != 2) { fprintf(stderr, "Usage: replay_features masked_outputs.csv\n"); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 2; }
    Campaign campaigns[32] = {0};
    int count = 0;
    char line[512];
    if (!fgets(line, sizeof(line), f)) return 2;
    while (fgets(line, sizeof(line), f)) {
        char name[96], group[16];
        int trace;
        double x, y;
        if (sscanf(line, "%95[^,],%d,%15[^,],%lf,%lf", name, &trace, group, &x, &y) != 5)
            return 2;
        if (strcmp(group, "fixed") && strcmp(group, "random")) return 2;
        int index;
        for (index = 0; index < count && strcmp(campaigns[index].name, name); index++);
        if (index == count) {
            if (count == 32) return 2;
            strcpy(campaigns[count++].name, name);
        }
        int g = !strcmp(group, "random");
        add(&campaigns[index].hw[g][0], x);
        add(&campaigns[index].hw[g][1], y);
        mixed_add(&campaigns[index].pair[g], x, y);
    }
    if (ferror(f)) return 2;
    fclose(f);
    puts("campaign,n_fixed,n_random,t_hw0,t_hw1,t_cp");
    for (int i = 0; i < count; i++) {
        Campaign *c = &campaigns[i];
        double t[2], mean[2], var[2];
        for (int s = 0; s < 2; s++) {
            Stats *a = &c->hw[0][s], *b = &c->hw[1][s];
            if (a->n < 2 || b->n < 2) return 2;
            t[s] = test(a->m, a->m2/(a->n-1), a->n, b->m, b->m2/(b->n-1), b->n);
        }
        Mixed *a = &c->pair[0], *b = &c->pair[1];
        double mx = (a->sx+b->sx)/(a->n+b->n), my = (a->sy+b->sy)/(a->n+b->n);
        mixed_z(a, mx, my, &mean[0], &var[0]);
        mixed_z(b, mx, my, &mean[1], &var[1]);
        printf("%s,%.0f,%.0f,%.17g,%.17g,%.17g\n", c->name, a->n, b->n, t[0], t[1],
               test(mean[0], var[0], a->n, mean[1], var[1], b->n));
    }
    return 0;
}
