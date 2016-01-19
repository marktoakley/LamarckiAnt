/*
    Parameters
    ----------
    paths : 2d array
        paths[t, r] is the position of replica r at time t

    Notes
    -----
    diffusion[p1, p2] : the average time elapsed since the replica at position
        p1 has been at p2.  time is measured in number of exchanges.

    rdiff[r, p] : the time elapsed since replica r was at position p



    for t, positions in enumerate(paths):
        maxelapsed += 1
        maxelapsed_cum += maxelapsed
        rdiff[:,:] += 1
        for r, p in enumerate(positions):
            rdiff[r, p] = 0

        for r in xrange(nreplicas):
            p1 = paths[t, r]
            for p2, telapsed in enumerate(rdiff[r,:]):
                diffusion[p1, p2] += telapsed
*/

#include <stdio.h>

int calculate_diffusion(long int * const paths, long int *diffusion, long int *rdiff,
    long int nreplicas, long int ntimes, long int *maxelapsed, long int *maxelapsed_cum)
{
  long int p, p1, p2, r, t, i, telapsed;
  long int *positions;

  for (i=0; i<nreplicas*nreplicas; ++i){
    rdiff[i] = 0;
    diffusion[i] = 0;
  }

  for (t=0; t<ntimes; ++t){
    *maxelapsed += 1;
    *maxelapsed_cum += *maxelapsed;

    for (i=0; i<nreplicas*nreplicas; ++i)
      rdiff[i] += 1;

    positions = &paths[t * nreplicas];
    for (r=0; r<nreplicas; ++r){
      p = positions[r];
      rdiff[r*nreplicas + p] = 0;
    }

    for (r=0; r<nreplicas; ++r){
      p1 = positions[r];
      for (p2=0; p2<nreplicas; ++p2){
        telapsed = rdiff[r*nreplicas + p2];
        diffusion[p1*nreplicas + p2] += telapsed;
      }
    }

  }
  return 0;
}
