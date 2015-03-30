/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package warpt;

/*
TODO: in preconditioner M, subtract horizontal average of shifts s
*/

import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.RecursiveExponentialFilter;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.Histogram;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.lapack.DMatrix;
import edu.mines.jtk.util.MedianFinder;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;


import dnp.CgSolver;
import dnp.InverseInterpolator;
import dnp.Vec;
import dnp.VecArrayFloat2;

/**
 * Dynamic warping for alignment of well logs. 
 * <p>
 * This application of dynamic warping differs from others in that it must
 * account for missing log values, and the fact that values in different logs
 * are measured at different depths.
 * <p>
 * The alignment of two log sequences f[i] and g[j] is represented by a
 * sequence of index pairs (i,j). This representation is the same as that
 * proposed by Sakoe and Chiba (1978) in their description of what is now
 * known as dynamic time warping (DTW). However, unlike Sakoe and Chiba, we
 * need not assume that the first and last samples of the sequence f[i] are
 * aligned with the first and last samples of the sequence g[j]. Indeed,
 * that assumption is rarely valid for well log sequences.
 * <p>
 * As for DTW, the first step is to compute alignment errors for all (i,j),
 * subject to the constraint that |j-i|&le;lmax. The difference l = j-i is
 * called lag (or shift), and one can use specified constraints on geologic
 * dip and distances between wells to compute the maximum lag lmax.
 * <p>
 * As noted above, conventional DTW assumes that lag l is zero for the first
 * and last samples of the optimal path. To permit the optimal path to begin
 * and end with non-zero lag, alignment errors are computed in a rotated
 * coordinate system: e[k,l] = pow(abs(f[i]-g[j]),epow), where k = j+i and l =
 * j-i. Here, k = imin+jmin,...,imax+jmax, and l = -lmax,...,lmax. (The
 * zero-based array index for any lag l is simply l+lmax.) When accumulating
 * alignment errors, we begin at kmin = imin+jmin and end at kmax = imax+jmax.
 * <p>
 * Half of the samples in the array of alignment errors e[k,l] are unused,
 * because i = (k-l)/2 and j = (k+l)/2 are integers for only even values of
 * k+l. For display purposes only, errors e[k,l] for which k+l is an odd
 * number may be computed by linear interpolation of the other errors.
 * <p>
 * Alignment errors e[k,l] for two log sequences f[i] and g[j] are computed
 * only for a range of k = i+j for which at least one of the two logs has a
 * non-null value. Within this range, where either f[i] or g[j] is null, the
 * null values are replaced by a non-null value randomly selected from the
 * sequence. This replacement does not actually alter the log sequences, as it
 * occurs only during the computation of alignment errors. Outside of this
 * range, alignment errors are set to a null error. When accumulating
 * alignment errors, any null errors are ignored. Accumulation of alignment
 * errors begins and ends with the first and last indices k for which
 * alignment errors are not null.
 *
 * @author Loralee Wheeler, Colorado School of Mines
 * @author Dave Hale, Colorado School of Mines
 * @version 2015.01.30
 */
public class WellLogWarping {

  /**
   * Sets the maximum shift (lag).
   * @param lmax the maximum lag.
   */
  public void setMaxShift(int lmax) {
    _lmax = lmax;
  }

  /**
   * Sets the exponent (the power) used to compute alignment errors.
   * @param epow the exponent.
   */
  public void setPowError(double epow) {
    _epow = (float)epow;
  }

  public void setPowError(double[] epow) {
     float[] ep = new float[epow.length];
     for (int i=0; i<epow.length; ++i) {
       ep[i] = (float)epow[i];
     }
    _epower = copy(ep);
  }

  /**
   * Sets the alignment error that represents no computed error. 
   * @param enull the null error.
   */
  public void setNullError(float enull) {
    _enull = enull;
  }

  /**
   * Sets the log value that represents no measured value.
   * The default null value is -999.2500.
   * @param vnull the null value.
   */
  public void setNullValue(float vnull) {
    _vnull = vnull;
  }

  /**
   * Returns an array of alignment errors e[k,l] for two sequences.
   * @param f array of values f[i] in 1st log sequence.
   * @param g array of values g[j] in 2nd log sequence.
   * @return array of alignment errors e[k,l].
   */
  public float[][] computeErrors(int m, float[] f, float[] g) {
    int ni = f.length;
    int nj = g.length;
    int lmax = min(max(ni,nj)-1,_lmax);
    int lmin = -lmax;
    int nl = 1+lmax-lmin;
    int nk = ni+nj-1;
    float[][] e = fillfloat(_enull,nl,nk);
    int[] igood = findGood(f);
    int[] jgood = findGood(g);
    int imin = min(igood);
    int imax = max(igood);
    int jmin = min(jgood);
    int jmax = max(jgood);
    int kmin = imin+jmin;
    int kmax = imax+jmax;
    //trace("computeErrors: kmin="+kmin+" kmax="+kmax);
    Random random = new Random(314159);
    for (int k=kmin; k<=kmax; ++k) {
      for (int l=lmin,ll=l-lmin; l<=lmax; ++l,++ll) {
        if ((k+l)%2==0) {
          int i = (k-l)/2;
          int j = (k+l)/2;
          float fi = value(random,igood,f,i);
          float gj = value(random,jgood,g,j);
          e[k][ll] = error(fi,gj,m);
        }
      }
    }
    return e;
  }

  /**
   * Returns accumulated errors d[k,l].
   * Any null errors in e[k,l] will be null in d[k,l].
   * @param e array of alignment errors e[k,l].
   * @return array of accumulated errors d[k,l].
   */
  public float[][] accumulateErrors(float[][] e) {
    int nk = e.length;
    int nl = e[0].length;
    int lmax = (nl-1)/2;
    int lmin = -lmax;
    float[][] d = fillfloat(_enull,nl,nk);

    // Range of k for which to accumulate errors.
    int kmin = kminNotNull(e);
    int kmax = kmaxNotNull(e);
    //trace("accumulateErrors: kmin="+kmin+" kmax="+kmax);

    // Special case: k = kmin.
    int k=kmin,km1,km2;
    for (int l=lmin,ll=l-lmin; l<=lmax; ++l,++ll) {
      if ((k+l)%2==0)
        d[k][ll] = e[k][ll];
    }

    // Special case: k = kmin+1.
    k = kmin+1; km1 = k-1;
    for (int l=lmin,ll=l-lmin,lm=ll-1,lp=ll+1; l<=lmax; ++l,++lm,++ll,++lp) {
      if ((k+l)%2==0) {
        float da = lm>=0?d[km1][lm]:FLT_MAX;
        float dc = lp<nl?d[km1][lp]:FLT_MAX;
        d[k][ll] = min(da,dc)+e[k][ll];
      }
    }

    // General case: k = kmin+2 to kmax.
    for (k=kmin+2,km1=k-1,km2=k-2; k<=kmax; ++k,++km1,++km2) {
      for (int l=lmin,ll=l-lmin,lm=ll-1,lp=ll+1; l<=lmax; ++l,++lm,++ll,++lp) {
        if ((k+l)%2==0) {
          float da = lm>=0?d[km1][lm]:FLT_MAX;
          float db =       d[km2][ll];
          float dc = lp<nl?d[km1][lp]:FLT_MAX;
          float dm = dmin(da,db,dc);
          d[k][ll] = dm+e[k][ll];
        }
      }
    }
    return d;
  }

  /**
   * Returns the optimal warping path as pairs of sample indices (k,l).
   * The pairs are returned as an array of two arrays, one array for the
   * indices k and the other array for the corresponding indices l.
   * @param d array of accumulated errors.
   * @return array[2][] containing indices (k,l). The array[0] will
   *  contain the indices k, and the array[1] will contain the indices l. The
   *  lengths of these two arrays will equal the number of pairs on the
   *  warping path; this number is unknown before the optimal path has been
   *  found.
   */
  public int[][] findWarping(float[][] d) {
    int nk = d.length;
    int nl = d[0].length;
    int lmax = (nl-1)/2;
    int lmin = -lmax;

    // Range of k for which to find the optimal warping path.
    int kmin = kminNotNull(d);
    int kmax = kmaxNotNull(d);
    //trace("findWarping: kmin="+kmin+" kmax="+kmax);

    // Initially empty arrays for (k,l) pairs.
    int nw = 0;
    int[] kw = new int[nk];
    int[] lw = new int[nk];

    // Find lag l with minimum accumulated error at k = kmax.
    int kp = kmax;
    int lp = -1;
    float dp = FLT_MAX;
    for (int l=lmin,ll=0; l<=lmax; ++l,++ll) {
      if ((kp+l)%2==0) {
        if (d[kp][ll]<dp) {
          lp = l;
          dp = d[kp][ll];
        }
      }
    }

    // Add the corresponding pair (k,l) to the path.
    kw[0] = kp;
    lw[0] = lp;
    nw += 1;

    // While the path is not yet complete, backtrack.
    while (kp>kmin) {
      int ll = lp-lmin;
      float da = lp>lmin  ?d[kp-1][ll-1]:FLT_MAX;
      float db = kp>kmin+1?d[kp-2][ll  ]:FLT_MAX;
      float dc = lp<lmax  ?d[kp-1][ll+1]:FLT_MAX;
      float dm = dmin(da,db,dc);
      if (dm==db) {
        kp -= 2;
      } else if (dm==da) {
        kp -= 1;
        lp -= 1;
      } else {
        kp -= 1;
        lp += 1;
      }
      kw[nw] = kp;
      lw[nw] = lp;
      nw += 1;
    }

    // Remove any wasted space from the arrays of indices, while reordering
    // the indices (k,l) so that k are increasing, not decreasing.
    int[] kt = new int[nw];
    int[] lt = new int[nw];
    for (int mw=0; mw<nw; ++mw) {
      kt[mw] = kw[nw-1-mw];
      lt[mw] = lw[nw-1-mw];
    }
    return new int[][]{kt,lt};
  }

  /**
   * Returns an optimal warping path converted from (k,l) to (i,j).
   * Omits any indices (i,j) for which either f[i] is null or g[j] is null.
   * @param kl array[2][] containing indices (k,l). The array[0] contains
   *  the indices k, and the array[1] contains the indices l.
   * @param f array of values f[i] in 1st log sequence.
   * @param g array of values g[j] in 2nd log sequence.
   * @return array[2][] containing indices (i,j). The array[0] contains
   *  the indices i, and the array[1] contains the indices j.
   */
  public int[][] convertWarping(int[][] kl, float[] f, float[] g) {
    int ni = f.length;
    int nj = g.length;
    int[] ks = kl[0];
    int[] ls = kl[1];
    int nkl = ks.length;

    // Initially empty arrays for index pairs (i,j).
    int nij = 0;
    int[] is = new int[nkl];
    int[] js = new int[nkl];

    // Collect index pairs (i,j) for all non-null values.
    for (int ikl=0; ikl<nkl; ++ikl) {
      int k = ks[ikl];
      int l = ls[ikl];
      int i = (k-l)/2;
      int j = (k+l)/2;
      //if (0<=i && i<ni && 0<=j && j<nj) {
      if (0<=i && i<ni && 0<=j && j<nj && f[i]!=_vnull && g[j]!=_vnull) {
        is[nij] = i;
        js[nij] = j;
        ++nij;
      }
    }

    // Return trimmed arrays of index pairs (i,j).
    is = copy(nij,is);
    js = copy(nij,js);
    return new int[][]{is,js};
  }

  /**
   * Apply warping to specified sequences f[i] and g[j].
   * The returned sequences f and g will both be indexed by k = i+j.
   * and will include null values for any missing values.
   * @param kl array[2][] containing indices (k,l). The array[0] contains
   *  the indices k, and the array[1] contains the indices l.
   * @param f array of values f[i] in 1st log sequence.
   * @param g array of values g[j] in 2nd log sequence.
   * @return array[2][] containing warped sequences f[k] and g[k].
   */
  public float[][] applyWarping(int[][] kl, float[] f, float[] g) {
    int[] ks = kl[0];
    int[] ls = kl[1];
    int nkl = ks.length;
    int ni = f.length;
    int nj = g.length;
    int nk = ni+nj-1;
    float[] fk = fillfloat(_vnull,nk);
    float[] gk = fillfloat(_vnull,nk);
    for (int ikl=0; ikl<nkl; ++ikl) {
      int k = ks[ikl];
      int l = ls[ikl];
      int i = (k-l)/2;
      int j = (k+l)/2;
      //fk[k] = (0<=i && i<ni && f[i]!=_vnull)?f[i]:_vnull;
      //gk[k] = (0<=j && j<nj && g[j]!=_vnull)?g[j]:_vnull;
      if (0<=i && i<ni && f[i]!=_vnull && 
          0<=j && j<nj && g[j]!=_vnull) {
        fk[k] = f[i];
        gk[k] = g[j];
      }
      if (ikl<nkl-1 && ks[ikl+1]==k+2) {
        fk[k+1] = fk[k];
        gk[k+1] = gk[k];
      }
    }
    return new float[][]{fk,gk};
  }

  /**
   * Interpolates alignment (or accumulated) errors for odd k+l.
   * Does not modify errors e[k,l] for which k+l is even.
   * <p>
   * Errors for odd k+l are never used. This interpolation is useful only for
   * displays of errors.
   * @param e array of errors e[k,l].
   */
  public void interpolateOddErrors(float[][] e) {
    int nk = e.length;
    int nl = e[0].length;
    int lmax = (nl-1)/2;
    int lmin = -lmax;
    for (int k=0; k<nk; ++k) {
      for (int l=lmin; l<=lmax; ++l) {
        if ((k+l)%2!=0) {
          int km = k-1; if (km<   0) km += 2;
          int kp = k+1; if (kp>= nk) kp -= 2;
          int lm = l-1; if (lm<lmin) lm += 2;
          int lp = l+1; if (lp>lmax) lp -= 2;
          float ekm = e[km][l-lmin];
          float ekp = e[kp][l-lmin];
          float elm = e[k][lm-lmin];
          float elp = e[k][lp-lmin];
          if (ekm==_enull || ekp==_enull || elm==_enull || elp==_enull) {
            e[k][l-lmin] = _enull;
          } else {
            e[k][l-lmin] = 0.25f*(ekm+ekp+elm+elp);
          }
        }
      }
    }
  }

  /**
   * Returns shifts for each specified log. Logs are assumed to have been
   * resampled so that every log has the same depth sampling. The returned
   * shifts are in units of samples, but may have non-zero fractional parts.
   * @param fs array[nl][nz] of log values, nz values for each of nl logs.
   * @return array[nl][nz] of shifts.
   */
  public float[][] findShifts(final float[] wm, final float[][][] wl) {
    int nz = wl[0][0].length;
    final int nm = wl[0].length;
    int nl = wl.length;
    int[] ls = new int[nl];
    final Pairs[] pt = new Pairs[nl*nm*(nl*nm-1)/2];
    final int[] ils = new int[nl*(nl-1)/2];
    final int[] jls = new int[nl*(nl-1)/2];
    int nlp = 0;
    for (int il=0; il<nl; ++il) {
      for (int jl=il+1; jl<nl; ++jl,++nlp) {
        ils[nlp] = il;
        jls[nlp] = jl;
      }
    }
    final AtomicInteger nij = new AtomicInteger(0);
    Parallel.loop(nlp,new Parallel.LoopInt() {
      public void compute(int ip) {
        for (int im=0; im<nm; ++im) {
          int il = ils[ip];
          int jl = jls[ip];
          float[] fi = wl[il][im];
          float[] gj = wl[jl][im];
          if (wm[im]>0 && wellNotNull(fi) && wellNotNull(gj)) {
            float[][] e = computeErrors(im,fi,gj);
            float[][] d = accumulateErrors(e);
            int[][] kl = findWarping(d);
            int[][] ij = convertWarping(kl,fi,gj);
            int[] is = ij[0];
            int[] js = ij[1];
            int np = is.length;
            if (np>0)
              pt[nij.getAndIncrement()] = new Pairs(il,jl,im,is,js,np);
          }
        }
      }
    });
    nlp = nij.get();
    Pairs[] ps = new Pairs[nlp];
    for (int ip=0; ip<nlp; ++ip)
      ps[ip] = pt[ip];
    computeWeights(wm,wl,ps,ls);

   return computeShifts(nz,ls,ps);
 }

  /**
   * Applies the specified shifts to the specified logs.
   * @param f array[nl][nz] of log values.
   * @param s array[nl][nz] of shifts.
   * @return array[nl][nz] of shifted log values.
   */
  public float[][] applyShifts(float[][] f, float[][] s) {
    int nk = f[0].length;
    int nl = f.length;
    float[][] g = fillfloat(_vnull,nk,nl);

    // For all logs, ...
    for (int il=0; il<nl; ++il) {

      // For all depths, ...
      for (int ik=0; ik<nk; ++ik) {

        // Depth (in samples) at which to interpolate log.
        //float zk = tts[ik]-s[il][ik];
        float zk = ik-s[il][ik];

        // Nearest-neighbor interpolation.
        int jk = (int)(zk+0.5f);
        if (0<=jk && jk<nk && f[il][jk]!=_vnull) {
          g[il][ik] = f[il][jk];
        }
      }
    }
    return g;
  }

  /**
   * Sorts well indices to approximately minimize distances between wells.
   * Specifically, this method approximately minimizes the total distance
   * traveled from well to well in a sequential iteration over all well
   * locations, in which each well is visited only once. Because exactly 
   * minimizing this total distance would require a costly solution to the
   * traveling-salesman problem, this method instead uses a simple greedy
   * solution.
   * <p>
   * This method is useful primarily in 2D displays of logs. Logs displayed as
   * adjacent pixels or curves are likely to appear more correlated than they
   * would be for an arbitrary ordering of wells.
   * @param x array of x coordinates of well locations.
   * @param y array of y coordinates of well locations.
   * @return the sorted array of well indices.
   */
  public static int[] sortWells(double[] x, double[] y) {
    int nw = x.length;
    Random r = new Random(314159);
    double dsmin = DBL_MAX; // the minimized distance
    int[] ksmin = null;
    for (int mw=0; mw<nw; ++mw) { // for all possible first-well indices
      boolean[] bs = new boolean[nw]; // flags for visited wells
      int[] ks = new int[nw]; // indices of visited wells
      int kw = mw; // index of the first well
      double xk = x[kw]; // x-coordinate
      double yk = y[kw]; // y-coordinate
      double ds = 0.0f; // distance sum
      int iw = 0;
      ks[iw] = kw; // first index in the list
      bs[kw] = true; // have now visited well with index kw
      for (iw=1; iw<nw; ++iw) { // for all other wells, ...
        int jmin = -1;
        double dmin = DBL_MAX;
        for (int jw=0; jw<nw; ++jw) { // find nearest not yet visited
          if (!bs[jw]) { // if not yet visited, ...
            double xj = x[jw];
            double yj = y[jw];
            double dx = xk-xj;
            double dy = yk-yj;
            double dj = dx*dx+dy*dy; // distance-squared
            if (dj<dmin) { // if nearest so far, ...
              dmin = dj;
              jmin = jw;
            }
          }
        }
        kw = jmin; // visit the nearest well
        xk = x[kw];
        yk = y[kw];
        ds += dmin;
        ks[iw] = kw;
        bs[kw] = true; // mark this well as visited
      }
      //trace("sortWells: ds="+ds);
      if (ds<dsmin) { // if this path has less distance, remember it
        dsmin = ds;
        ksmin = ks;
      }
    }
    //trace("sortWells: dsmin="+dsmin);
    return ksmin;
  }

  /**
   * Gets a uniform sampling of depth for specified depths and logs.
   * Considers only depths for which log values are non-null. The returned
   * sampling will include the shallowest and the deepest depths logged.
   * @param z array of arrays of depths; one array for each log.
   * @param f array of arrays of log values; one array for each log.
   * @return the uniform sampling.
   */
  public Sampling getDepthSampling(float[][] z, float[][] f) {
    int nl = z.length;

    // Total number of depths specified.
    int nlz = 0;
    for (int il=0; il<nl; ++il)
      nlz += z[il].length;

    // Array for depth increments, and counter for number of increments.
    float[] dz = new float[nlz];
    int ndz = 0;

    // Find min and max depths, while storing depth increments.
    // Consider only log samples with non-null values.
    float zmin =  FLT_MAX;
    float zmax = -FLT_MAX;
    for (int il=0; il<nl; ++il) {
      int nz = z[il].length;
      float zi = z[il][0];
      float fi = f[il][0];
      for (int iz=0; iz<nz; ++iz) {
        float zim1 = zi;
        float fim1 = fi;
        zi = z[il][iz];
        fi = f[il][iz];
        if (fi!=_vnull && zi<zmin)
          zmin = zi;
        if (fi!=_vnull && zi>zmax)
          zmax = zi;
        if (iz>0 && fi!=_vnull && fim1!=_vnull)
          dz[ndz++] = zi-zim1;
      }
    }

    // Depth interval is median of all depth increments.
    dz = copy(ndz,dz);
    MedianFinder mf = new MedianFinder(ndz);
    float zdel = mf.findMedian(dz);

    // Uniform sampling.
    int nz = 1+(int)ceil((zmax-zmin)/zdel);
    return new Sampling(nz,zdel,zmin);
  }

  public Sampling getDepthSampling(float[][][] z, float[][][] f) {
    int nw = z.length;
    int nl = z[0].length;

    // Total number of depths specified.
    int nlz = 0;

    for (int iw=0; iw<nw; ++iw)
      for (int il=0; il<nl; ++il)
        nlz += z[iw][il].length;

    // Array for depth increments, and counter for number of increments.
    float[] dz = new float[nlz];
    int ndz = 0;

    // Find min and max depths, while storing depth increments.
    // Consider only log samples with non-null values.
    float zmin =  FLT_MAX;
    float zmax = -FLT_MAX;
    for (int iw=0; iw<nw; ++iw) {
      for (int il=0; il<nl; ++il) {
        int nz = z[iw][il].length;
        if (nz!=0) {
          float zi = z[iw][il][0];
          float fi = f[iw][il][0];
          for (int iz=0; iz<nz; ++iz) {
            float zim1 = zi;
            float fim1 = fi;
            zi = z[iw][il][iz];
            fi = f[iw][il][iz];
            if (fi!=_vnull && zi<zmin)
              zmin = zi;
            if (fi!=_vnull && zi>zmax)
              zmax = zi;
            if (iz>0 && fi!=_vnull && fim1!=_vnull)
              dz[ndz++] = zi-zim1;
          }
        }
      }
    }

    // Depth interval is median of all depth increments.
    dz = copy(ndz,dz);
    MedianFinder mf = new MedianFinder(ndz);
    float zdel = mf.findMedian(dz);

    // Uniform sampling.
    int nz = 1+(int)ceil((zmax-zmin)/zdel);
    return new Sampling(nz,zdel,zmin);
  }

  /**
   * Resamples a log with specified depths and values.
   * The desired depth sampling for the output array of values can be
   * different from (e.g., coarser than) that implied by the input array of
   * depths.
   * @param s the desired uniform sampling of depths.
   * @param z array of depths for which values are provided.
   * @param f array of values; some values may be null.
   * @return array of values for uniformly sampled depths.
   */
  public float[] resampleLog(Sampling s, float[] z, float[] f) {
    int nz = s.getCount();
    float zmin = (float)s.getFirst();
    float zmax = (float)s.getLast();
    int n = z.length;
    float[] g = new float[nz];
    float[] c = new float[nz];
    if (n!=0) {
      for (int i=0; i<n; ++i) {
        float zi = z[i];
        float fi = f[i];
        if (zmin<=zi && zi<=zmax && fi!=_vnull) {
          int iz = s.indexOfNearest(zi);
          g[iz] += fi;
          c[iz] += 1.0f;
        }
      }
      for (int iz=0; iz<nz; ++iz) {
        if (c[iz]>0.0f) {
          g[iz] /= c[iz];
        } else {
          g[iz] = _vnull;
        }
      }
    } else {
      g = fillfloat(_vnull,nz);
    }
    return g;
  }

  /**
   * Resamples multiple logs with specified depths and values.
   * This method simply calls the method 
   * {@link #resampleLog(Sampling,float[],float[])}
   * for all specified arrays of depths and values.
   * @param s the desired uniform sampling of depths.
   * @param z array of depths for which values are provided.
   * @param f array of values; some values may be null.
   * @return array of values for uniformly sampled depths.
   */
  public float[][] resampleLogs(Sampling s, float[][] z, float[][] f) {
    int n = z.length;
    float[][] g = new float[n][];
    for (int i=0; i<n; ++i)
      g[i] = resampleLog(s,z[i],f[i]);
    return g;
  }
  public float[][][] resampleLogs(Sampling s, float[][][] z, float[][][] f) {
    int n = z.length;
    float[][][] g = new float[n][][];
    for (int i=0; i<n; ++i)
      g[i] = resampleLogs(s,z[i],f[i]);
    return g;
  }
  public void replaceNullsS(float[] f, float freplace, float[] s) {
    int n = f.length;
    for (int i=0; i<n; ++i) {
      if (f[i]==_vnull) {
        s[i] = freplace;
      } 
    }
  }
  public void replaceNullsS(float[][] f, float freplace, float[][] s) {
    int n = f.length;
    float[][] g = new float[n][];
    for (int i=0; i<n; ++i)
      replaceNullsS(f[i],freplace,s[i]);
  }

  /**
   * Replaces any null values found in a log with a specified value.
   * Typically used only for displays.
   * @param f array of log values.
   * @param freplace value used to replace any null values.
   * @return array of logs values with null values replaced.
   */
  public float[] replaceNulls(float[] f, float freplace) {
    int n = f.length;
    float[] g = new float[n];
    for (int i=0; i<n; ++i) {
      if (f[i]!=_vnull) {
        g[i] = f[i];
      } else {
        g[i] = freplace;
      }
    }
    return g;
  }
  /**
   * Replaces any null values found in multiple logs with a specified value.
   * Typically used only for displays.
   * @param f array of log values.
   * @param freplace value used to replace any null values.
   * @return array of logs values with null values replaced.
   */
  public float[][] replaceNulls(float[][] f, float freplace) {
    int n = f.length;
    float[][] g = new float[n][];
    for (int i=0; i<n; ++i)
      g[i] = replaceNulls(f[i],freplace);
    return g;
  }

  public float[] toFloat(int[] i) {
    int n = i.length;
    float[] f = new float[n];
    for (int j=0; j<n; ++j)
      f[j] = (float)i[j];
    return f;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
 
  private int _lmax = Integer.MAX_VALUE;
  private float _epow = 1.0f;
  private float[] _epower = null;
  private float _enull = -FLT_MIN;
  private float _vnull = -999.2500f;
  private static final float SUM_SCL = 0.001f;

  /**
   * Arrays of pairs of depth sample indices (is,js) with weights ws. These
   * arrays were computed by warping a pair of logs with indices (ilog,jlog).
   */
  private static class Pairs {
    Pairs(int il, int jl, int m, int[] is, int[] js, int np)
  {
      this.il = il;
      this.jl = jl;
      this.is = is;
      this.js = js;
      this.np = np;
      this.m  = m;
    }
    int il,jl,m,np;
    int[] is,js;
    double ws;
  }

  /**
   * Conjugate-gradient operator A and preconditioner M.
   * The preconditioner smooths along depth, while subtracting
   * the mean and linear trend.
   */
  private static class A implements CgSolver.A {
    A(float[][] tz, int[][] c, double wp, Pairs[] ps) {
      _ps = ps;
      _tz = tz;
      _tc = c;
      _wp = wp;
    }
    public void apply(Vec vx, Vec vy) {
      float[][] x = ((VecArrayFloat2)vx).getArray();
      float[][] y = ((VecArrayFloat2)vy).getArray();
      applyLhs(_tz,_tc,_wp,_ps,x,y);
      testNND(y[0].length,y.length);
    }
    private void testNND(int nk, int nl) {
      float[][] x  = sub(randfloat(nk,nl),0.5f);
      float[][] xt = copy(x);
      applyLhs(_tz,_tc,_wp,_ps,copy(x),x);

      VecArrayFloat2 vx = new VecArrayFloat2(x);
      VecArrayFloat2 vxt = new VecArrayFloat2(xt);

      float[][] y  = sub(randfloat(nk,nl),0.5f);
      float[][] yt = copy(y);
      applyLhs(_tz,_tc,_wp,_ps,copy(y),y);

      VecArrayFloat2 vy  = new VecArrayFloat2(y);
      VecArrayFloat2 vyt = new VecArrayFloat2(yt);

      float ytAx = (float)vyt.dot(vx);
      float xtAy = (float)vxt.dot(vy);
      float xtAx = (float)vxt.dot(vx);
      float ytAy = (float)vyt.dot(vy);
      //System.out.println("ytAx="+ytAx+" xtAy="+xtAy);
      //assert abs(ytAx-xtAy)<1.0e-8 : "A is not symmetric";
      assert (xtAx>=0 || ytAy>=0) : "A is not non-negative";
    }
    private Pairs[] _ps;
    private float[][] _tz;
    private int[][] _tc;
    private double _wp;
  }
  private static class Ml implements CgSolver.A {
    Ml(double sigma, int nk, int nl, int[] ls) {
      _ref = new RecursiveExponentialFilter(sigma);
      _ref.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE);
      _s = new float[nk];
      _ls = ls;
    }
    public void apply(Vec vx, Vec vy) {
      float[][] x = ((VecArrayFloat2)vx).getArray();
      float[][] y = ((VecArrayFloat2)vy).getArray();
      copy(x,y);
      subtractMeanOverLogs(y);
      _ref.apply1(y,y);
      subtractMeanOverLogs(y);
    }
    public void subtractMeanOverLogs(float[][] x) {
      int nk = x[0].length;
      int nl = x.length;
      int cl = 0;
      zero(_s);
      for (int il=0; il<nl; ++il) {
        if (_ls[il]>0) {
          for (int ik=0; ik<nk; ++ik)
            _s[ik] += x[il][ik];
          ++cl;
        }
      }
      float c = 1.0f/cl;
      for (int il=0; il<nl; ++il) {
        if (_ls[il]>0) {
          for (int ik=0; ik<nk; ++ik)
            x[il][ik] -= _s[ik]*c;
        }
      }
    }
    public void test(float[][] x) {
      int nk = x[0].length;
      int nl = x.length;
      zero(_s);
      for (int il=0; il<nl; ++il)
        for (int ik=0; ik<nk; ++ik)
          if (_ls[il]>0)
            _s[ik] += x[il][ik];
      trace("M.test: sum="); dump(_s);
    }
    private RecursiveExponentialFilter _ref;
    float[] _s; // used to efficiently compute sum over logs
    int[][] _tc;
    int[] _ls;
  }

  private static void makeRhs(
    float[][] tz, int[][] c, float[] sc, Pairs[] ps, float[][] y) 
  {
    int nt = y[0].length; // number of depths
    int nl = y.length; // number of logs
    int np = ps.length; // number of log pairs
    zero(y); // zero y before accumulating below
    zero(c);
    for (int ip=0; ip<np; ++ip) { // for all log pairs, ...
      Pairs p = ps[ip];
      int il = p.il;
      int jl = p.jl;
      float[] yi = y[il];
      float[] yj = y[jl];
      float[] ti = tz[il];
      float[] tj = tz[jl];
      int[] ci = c[il];
      int[] cj = c[jl];
      int[] is = p.is;
      int[] js = p.js;
      double ws = p.ws;
      int n = is.length;
      for (int k=0; k<n; ++k) { // for all index pairs (i,j), ...
        int ik = is[k];
        int jk = js[k];
        int it = (int)(ti[ik]+0.5f);
        int jt = (int)(tj[jk]+0.5f);
        if (0<it && it<nt && 0<jt && jt<nt) {
          float scl = (float)(ws*ws);
          float dif = jk-ik;
          ++ci[it];
          ++cj[jt];
          dif *= scl;
          yi[it] += dif;
          yj[jt] -= dif;
        }
      }
    }
  }

  private static void applyLhs(
    float[][] tz, int[][] c, double wp, Pairs[] ps, float[][] x, float[][] y) 
  {
    int nt = x[0].length; // number of depths
    int nl = x.length; // number of logs
    int np = ps.length; // number of log pairs
    zero(y); // zero y before accumulating below
    for (int ip=0; ip<np; ++ip) { // for all log pairs, ...
      Pairs p = ps[ip];
      int il = p.il;
      int jl = p.jl;
      float[] xi = x[il];
      float[] xj = x[jl];
      float[] yi = y[il];
      float[] yj = y[jl];
      float[] ti = tz[il];
      float[] tj = tz[jl];
      int[] is = p.is;
      int[] js = p.js;
      double ws = p.ws;
      int n = is.length;
      for (int k=0; k<n; ++k) { // for all index pairs (i,j), ...
        int ik = is[k];
        int jk = js[k];
        int it = (int)(ti[ik]+0.5f);
        int jt = (int)(tj[jk]+0.5f);
        if (0<it && it<nt && 0<jt && jt<nt) {
          float scl = (float)(ws*ws);
          float dif = xi[it]-xj[jt];
          dif *= scl;
          yi[it] += dif;
          yj[jt] -= dif;
        }
      }
    }
    constrainSmoothShifts(wp,c,x,y);
  }

  private static float dmin(float a, float b, float c) {
    float d = b;
    if (a<d) d = a;
    if (c<d) d = c;
    return d;
  }

  private float error(float f, float g, int m) {
    float d = f-g;
    float epow = _epower[m];
    if (epow==2.0f) {
      return d*d;
    } else if (epow==1.0f) {
      return abs(d);
    } else {
      if (d<0.0f) d = -d;
      return pow(d,epow);
    }
  }
  private float error(float f, float g) {
    float d = f-g;
    if (_epow==2.0f) {
      return d*d;
    } else if (_epow==1.0f) {
      return abs(d);
    } else {
      if (d<0.0f) d = -d;
      return pow(d,_epow);
    }
  }

  public int[] findGood(float[] f) {
    int n = f.length;
    int[] igood = new int[n];
    int ngood = 0;
    for (int i=0; i<n; ++i) {
      if (f[i]!=_vnull) {
        igood[ngood] = i;
        ++ngood;
      }
    }
    return copy(ngood,igood);
  }

  public float value(Random random, int[] igood, float[] f, int i) {
    int n = f.length;
    if (0<=i && i<n && f[i]!=_vnull) {
      return f[i];
    } else {
      int j = random.nextInt(igood.length);
      return f[igood[j]];
    }
  }

  private int kminNotNull(float[][] e) {
    int nk = e.length;
    int nl = e[0].length;
    for (int k=0; k<nk; ++k) {
      if (e[k][0]!=_enull || e[k][1]!=_enull)
        return k;
    }
    return nk;
  }

  private int kmaxNotNull(float[][] e) {
    int nk = e.length;
    int nl = e[0].length;
    for (int k=nk-1; k>=0; --k) {
      if (e[k][0]!=_enull || e[k][1]!=_enull)
        return k;
    }
    return -1;
  }

  // Post-processing and inversion of computed shifts.
  private static void cleanShifts(float[] s) {
    int n1 = s.length;
    for (int i1=1; i1<n1; ++i1) {
      if (s[i1]<s[i1-1]-0.999f)
        s[i1] = s[i1-1]-0.999f;
    }
  }
  private static void cleanShiftsR(float[] s) {
    int n1 = s.length;
    for (int i1=1; i1<n1; ++i1) {
      if (s[i1]>s[i1-1]+0.999f)
        s[i1] = s[i1-1]+0.999f;
    }
  }
  private static void invertShiftsR(float[] u, float[] t, float[] s) {
    cleanShiftsR(s);
    int n1 = s.length;
    for (int i1=0; i1<n1; ++i1)
      s[i1] = u[i1] - s[i1];
    inverseInterpolate(u,s,t);
    for (int i1=0; i1<n1; ++i1) 
      s[i1] = t[i1]-u[i1];
  }
  private static void invertShiftsR(float[][] s) {
    int n1 = s[0].length;
    int n2 = s.length;
    float[] u = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    for (int i2=0; i2<n2; ++i2)
      invertShiftsR(u,t,s[i2]);
  }
  private static void invertShifts(float[] u, float[] t, float[] s) {
    cleanShifts(s);
    int n1 = s.length;
    for (int i1=0; i1<n1; ++i1)
      s[i1] += u[i1];
    inverseInterpolate(u,s,t);
    for (int i1=0; i1<n1; ++i1) 
      s[i1] = u[i1]-t[i1];
  }
  private static void invertShifts(float[][] s) {
    int n1 = s[0].length;
    int n2 = s.length;
    float[] u = rampfloat(0.0f,1.0f,n1);
    float[] t = zerofloat(n1);
    for (int i2=0; i2<n2; ++i2)
      invertShifts(u,t,s[i2]);
  }
  public static void stats(float[] wl) {
    int nk = wl.length;
    int nn = 0;
    for (int k=0; k<nk; ++k) {
      if (wl[k]!=-999.2500) {
        ++nn;
      }
    }
    if (nn>0) {
      float[] good = zerofloat(nn);
      nn = 0;
      for (int k=0; k<nk; ++k) {
        if (wl[k]!=-999.2500) {
          good[nn] = wl[k];
          ++nn;
        }
      }
      int p25 = (int)(ceil(0.25*nn));
      int p50 = (int)(ceil(0.50*nn));
      int p75 = (int)(ceil(0.75*nn));
      quickPartialSort(p25,good);
      float g25 = good[p25];
      quickPartialSort(p50,good);
      float med = good[p50];
      quickPartialSort(p75,good);
      float g75 = good[p75];
      float iqr = (g75-g25);
      float mean = sum(good)/nn;
      float[] diff = sub(good,mean);
      float[] sq = mul(diff,diff);
      float stdev = sqrt(sum(sq)/nn);
      System.out.println("mean="+mean+" stdev="+stdev);
      System.out.println(" med="+med+"  iqr="+iqr);
    } else {
      System.out.println("Null log");
    }
  }

  private void clipShifts(float[][] s) {
    int nl = s.length;
    int nk = s[0].length;
    for (int l=0; l<nl; ++l) {
      for (int k=0; k<nk; ++k) {
        if (s[l][k]>_lmax) s[l][k] = _lmax;
        if (s[l][k]<-_lmax) s[l][k] = -_lmax;
      }
    }
  }
  
  private static void interpolate(float[] u, float[] x, float[] v, float[] y) {
    CubicInterpolator ci =
      new CubicInterpolator(CubicInterpolator.Method.LINEAR,u,x);
    ci.interpolate(v,y);
  }

  private static void inverseInterpolate(float[] u, float[] x, float[] y) {
    CubicInterpolator ci =
      new CubicInterpolator(CubicInterpolator.Method.LINEAR,x,u);
    ci.interpolate(u,y);
  }

  public boolean wellNotNull(float[] f) {
    int nk = f.length;
    for (int ik=0; ik<nk; ++ik) {
      if (f[ik]!=_vnull) return true;
    }
    return false;
  }

  private static float[] findConstantShift(int[] ls, Pairs[] ps) {
    int nl = ls.length;
    int nlp = ps.length; // number of log pairs
    int np = 1;
    for (int ip=0; ip<nlp; ++ip) // for all log pairs, ...
      np += ps[ip].np;
    DMatrix a = new DMatrix(np,nl);
    DMatrix b = new DMatrix(np,1);
    for (int il=0; il<nl; ++il) {
      if (ls[il]>0)
        a.set(0,il,1.0d); // makes sum = 0
    }
    np = 1;
    for (int ip=0; ip<nlp; ++ip) { // for all log pairs, ...
      Pairs p = ps[ip];
      int il = p.il;
      int jl = p.jl;
      int[] is = p.is;
      int[] js = p.js;
      double ws = p.ws;
      int n = is.length;
      for (int k=0; k<n; ++k,++np) { 
        int ik = is[k];
        int jk = js[k];
        double w = ws;
        double dz = (double)(jk-ik);
        a.set(np,il, ws);
        a.set(np,jl,-ws);
        b.set(np, 0, w*dz);
      }
    }
    int nls = 0;
    for (int il=0; il<nl; ++il)
      if (ls[il]>0)
        ++nls;
    DMatrix at = new DMatrix(np,nls);
    nls = 0;
    for (int il=0; il<nl; ++il) {
      if (ls[il]>0) {
        for (int ip=0; ip<np; ++ip)
          at.set(ip,nls,a.get(ip,il));
        ++nls;
      }
    }

    DMatrix x = at.solve(b);
    double[] sd = x.getArray();
    float[] s  = new float[nl];
    nls = 0;
    for (int il=0; il<nl; ++il) {
      if (ls[il]>0) {
        s[il] = (float)sd[nls];
        ++nls;
      }
    }

    return s;
  }

  // Use CG to solve least-squares equations for shifts s.
  private float[][] computeShifts(int nt, int[] ls, Pairs[] ps) {
   int nl = ls.length;
   int niter = 5;
   float[][] r = new float[nl][nt]; // r = s+q
   float[][] q = new float[nl][nt]; // handles shifts for str/sqz
   float[][] tz = new float[nl][nt]; 
   float[] s = findConstantShift(ls,ps); // static shifts
   //trace("sum="+sum(s));
   //trace("Static shifts logs:");
   //dump(s);
   for (int il=0; il<nl; ++il) { 
     tz[il] = rampfloat(0.0f,1.0f,nt); // initialize t = z + s
     add(tz[il],s[il],tz[il]);
     r[il] = fillfloat(s[il],nt);
   }
   
   for (int iter=0; iter<niter; ++iter) {
     innerLoop(iter,s,ls,ps,tz,q);
     trace("outer it="+iter);

     for (int il=0; il<nl; ++il) 
       add(q[il],s[il],r[il]);
     clipShifts(r);
     updateTz(tz,r);
     subtractMeanOverLogs(ls,r);
   }
   return r;
  }

  private void innerLoop(
    int iter, float[] s, int[] ls, Pairs[] ps, float[][] tz, float[][] q) 
  {
   int nt = tz[0].length;
   int nl = tz.length;
   int[][] c = new int[nl][nt]; // for count of used pairs

   float sigma = 100.0f;
   float small = 1.0e-6f;
   int niter = 200*iter+200;
   double wp = sumWeights(tz,ps);
   A a = new A(tz,c,wp,ps);
   Ml m = new Ml(sigma,nt,nl,ls);
   CgSolver cs = new CgSolver(small,niter);
   float[][] b = new float[nl][nt]; // for right-hand side
   float[][] p = new float[nl][nt]; 
   for (int il=0; il<nl; ++il)
     p[il] = fillfloat(s[il],nt);

   makeRhs(tz,c,s,ps,b);
   applyLhs(tz,c,wp,ps,copy(p),p);
   sub(b,p,b);

   VecArrayFloat2 vb = new VecArrayFloat2(b);
   VecArrayFloat2 vq = new VecArrayFloat2(q);
   //cs.solve(a,vb,vq); // gives q
   cs.solve(a,m,vb,vq); // gives q
   //m.test(q);
  }

  private double sumWeights(float[][] tz, Pairs[] ps) {
    int nt = tz[0].length;
    int np = ps.length; // number of log pairs
    double wp = 0;
    for (int ip=0; ip<np; ++ip) { // for all log pairs, ...
      Pairs p = ps[ip];
      int n = p.np;
      int il = p.il;
      int jl = p.jl;
      int[] is = p.is;
      int[] js = p.js;
      float[] ti = tz[il];
      float[] tj = tz[jl];
      double ws = p.ws;
      for (int k=0; k<n; ++k) { // for all index pairs (i,j), ...
        int ik = is[k];
        int jk = js[k];
        int it = (int)(ti[ik]+0.5f);
        int jt = (int)(tj[jk]+0.5f);
        if (0<it && it<nt && 0<jt && jt<nt)
          wp += ws;
      }
    }
    return wp;
  }

  private void updateTz(float[][] tz, float[][] r) {
    int nt = tz[0].length;
    int nl = tz.length;
    float[] ut = rampfloat(0.0f,1.0f,nt);
    float[] uz = rampfloat(0.0f,1.0f,nt);
    float[] zt = new float[nt];
    for (int il=0; il<nl; ++il) {
      cleanShiftsR(r[il]);
      for (int it=0; it<nt; ++it) 
        zt[it] = it-r[il][it];
      interpolate(zt,ut,uz,tz[il]); // inverts z(t) to t(z)
    }
  }

  private static void constrainSmoothShifts(
    double wp, int[][] c, float[][] x, float[][] y) 
  {
    int nt = x[0].length; // number of depths
    int nl = x.length; // number of logs
    wp = 400;
    for (int il=0; il<nl; ++il) {
      float[] xi = x[il];
      float[] yi = y[il];
      for (int it=1; it<nt; ++it) { 
        double wk = (c[il][it]==0 || c[il][it-1]==0) ? wp : 0.0d;
        float scl = (float)(wk*wk);
        float dif = xi[it]-xi[it-1];
        dif *= scl;
        yi[it  ] += dif;
        yi[it-1] -= dif;
      }
    }
  }

  public float medianAbsoluteDeviation(float[][] f) {
    int nt = f[0].length;
    int nl = f.length;
    float md;
    float[] mad = new float[nt*nl];
    int im = 0;

    for (int it=0; it<nt; ++it) {
      float[] med = new float[nl];
      int jl = 0;
      for (int il=0; il<nl; ++il) {
        if (f[il][it]!=_vnull)
          med[jl++] = f[il][it];
      }
      if (jl>0) {
        float[] medt = copy(jl,med);
        Arrays.sort(medt);
        if (jl%2==0)
          md = (medt[jl/2] + medt[jl/2-1])/2.0f;
        else
          md = medt[jl/2];
        jl = 0;
        for (int il=0; il<nl; ++il) {
          if (f[il][it]!=_vnull) {
            mad[im++] = abs(f[il][it]-md);
          }
        }
      }
    }
    float[] mad2 = copy(im,mad);
    Arrays.sort(mad2);
    if (im%2==0)
      md = (mad2[im/2] + mad2[im/2-1])/2.0f;
    else
      md = mad2[im/2];

    return md;
  }

  public float[][] replaceWithMAD(float[][] f) {
    int nt = f[0].length;
    int nl = f.length;
    float md;
    float[][] fm = fillfloat(_vnull,nt,nl);

    for (int it=0; it<nt; ++it) {
      float[] med = new float[nl];
      int jl = 0;
      for (int il=0; il<nl; ++il) {
        if (f[il][it]!=_vnull)
          med[jl++] = f[il][it];
      }
      if (jl>0) {
        float[] medt = copy(jl,med);
        Arrays.sort(medt);
        if (jl%2==0)
          md = (medt[jl/2] + medt[jl/2-1])/2.0f;
        else
          md = medt[jl/2];
        jl = 0;
        for (int il=0; il<nl; ++il) {
          if (f[il][it]!=_vnull) {
            fm[il][it] = abs(f[il][it]-md);
          }
        }
      }
    }
    return fm;
  }

  public float[][] replaceWithMedian(float[][] f) {
    int nt = f[0].length;
    int nl = f.length;
    float md;
    float[][] fm = fillfloat(_vnull,nt,nl);

    for (int it=0; it<nt; ++it) {
      float[] med = new float[nl];
      int jl = 0;
      for (int il=0; il<nl; ++il) {
        if (f[il][it]!=_vnull)
          med[jl++] = f[il][it];
      }
      if (jl>0) {
        float[] medt = copy(jl,med);
        Arrays.sort(medt);
        if (jl%2==0) {
          md = (medt[jl/2] + medt[jl/2-1])/2.0f;
        } else {
          md = medt[jl/2];
        }
        for (int il=0; il<nl; ++il) {
          if (f[il][it]!=_vnull)
            fm[il][it] = md;
        }
      }
    }
    return fm;
  }

  private void computeWeights(float[] wm, float[][][] wl, Pairs[] ps, int[] ls) {
    int  nm = wl[0].length;
    int nlp = ps.length;
    double[] wsum = new double[nm];
    for (int lp=0; lp<nlp; ++lp) {
      int  m = ps[lp].m;
      int np = ps[lp].np;
      int il = ps[lp].il;
      int jl = ps[lp].jl;
      int[] is = ps[lp].is;
      int[] js = ps[lp].js;
      double ws = ps[lp].ws;
      float tae = 0;
      for (int kz=0; kz<np; ++kz) 
        tae += error(wl[il][m][is[kz]],wl[jl][m][js[kz]],m);
      ws = tae>0 ? np*pow(np/tae,2.0f/_epower[m]) : 1.0f; 
      ps[lp].ws = ws;
      wsum[m] += ws;
      ls[il] += 1;
      ls[jl] += 1;
    }
    for (int lp=0; lp<nlp; ++lp) {
      int m = ps[lp].m;
      double ws = ps[lp].ws;
      ws /= wsum[m];
      ws *= wm[m];
      //ps[lp].ws = 1.0f;
      ps[lp].ws = ws;
      trace("il="+ps[lp].il+" jl="+ps[lp].jl+" w="+ws);
    }
  }


  private static void trace(String s) {
    System.out.println(s);
  }
  private void subtractMeanOverLogs(int[] ls, float[][] x) {
    int nk = x[0].length;
    int nl = x.length;
    int cl = 0;
    float[] s = zerofloat(nk);
    for (int il=0; il<nl; ++il) {
      if (ls[il]>0) {
        for (int ik=0; ik<nk; ++ik)
          s[ik] += x[il][ik];
        ++cl;
      }
    }
    float c = 1.0f/cl;
    for (int il=0; il<nl; ++il) {
      if (ls[il]>0) {
        for (int ik=0; ik<nk; ++ik)
          x[il][ik] -= s[ik]*c;
      }
    }
  }
}
