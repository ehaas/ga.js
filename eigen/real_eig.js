(function() {
  var util = require('util');
  function sqmrOut(N, matrix, sPar){

  if (N > 0) {

    sPar.outmsgstr = sPar.outmsgstr + "<br /> &nbsp; &nbsp;";

    for (var i = 0; i < N; i++) {
      for (var j = 0; j < N; j++)  sPar.outmsgstr = sPar.outmsgstr + matrix[i][j] + " &nbsp; &nbsp; ";
      sPar.outmsgstr = sPar.outmsgstr + "<br /> &nbsp; &nbsp;";
    } // End for i

  }   // End if (N > 0)

    return;
  }  //End of sqmrOut


  function outputcomplexvector(erFlag, vecSize, rVec, iVec, sPar){

  if (vecSize > 0) {

    sPar.outmsgstr = sPar.outmsgstr + "<br /> &nbsp; &nbsp;";

    for (var i = 0; i < vecSize; i++) {
      sPar.outmsgstr = sPar.outmsgstr + rVec[i];
      if (iVec[i] != 0) {

        if (iVec[i] > 0) {
          sPar.outmsgstr = sPar.outmsgstr + " &nbsp; + &nbsp; " + iVec[i] + " i" ;
        }  // End if (iVec[i] > 0)
        else {
          sPar.outmsgstr = sPar.outmsgstr + " &nbsp; - &nbsp; " + Math.abs(iVec[i]) + " i" ;
        }  // End else (iVec[i] < 0)
        
      }  // End if (iVec[i] != 0)
      
      sPar.outmsgstr = sPar.outmsgstr + "<br /> &nbsp; &nbsp;";

    } // End for i
    
    if (erFlag != 99){ // If erFlag equals 99, do not print Error Code
      sPar.outmsgstr = sPar.outmsgstr + "<br /> &nbsp; &nbsp; <b>Error Code</b>: " + erFlag;
    } // End if (erFlag != 99)

    sPar.outmsgstr = sPar.outmsgstr + "<br />";
  }   // End if (vecSize > 0)

    return;
  }  //End of outputcomplexvector


  function cdivA(ar, ai, br, bi, A, in1, in2, in3){
  // Division routine for dividing one complex number into another:
  // This routine does (ar + ai)/(br + bi) and returns the results in the specified
  // elements of the A matrix.

   var s, ars, ais, brs, bis;

   s = Math.abs(br) + Math.abs(bi);
   ars = ar/s;
   ais = ai/s;
   brs = br/s;
   bis = bi/s;
   s = brs*brs + bis*bis;
   A[in1][in2] = (ars*brs + ais*bis)/s;
   A[in1][in3] = (-(ars*bis) + ais*brs)/s;
   return;
  } // End cdivA

  function hqr2(N, A, B, low, igh, wi, wr, oPar){
  /* Computes the eigenvalues and eigenvectors of a real upper-Hessenberg Matrix using the QR method. */

  var norm = 0.0, p, q, ra, s, sa, t = 0.0, tst1, tst2, vi, vr, w, x, y, zz;
  var k = 0, l, m, mp2, en = igh, incrFlag = 1, its, itn = 30*N, enm2, na, notlas;

  for (var i = 0; i < N; i++){ // Store eigenvalues already isolated and compute matrix norm.
   for (var j = k; j < N; j++)
    norm += Math.abs(A[i][j]);
   k = i;
   if ((i < low) || (i > igh)){
    wi[i] = 0.0;
    wr[i] = A[i][i];
   } //End if (i < low or i > igh)
  }  // End for i

  // Search next eigenvalues
  while (en >= low){

   if (incrFlag) { //Skip this part if incrFlag is set to 0 at very end of while loop
    its = 0;
    na = en - 1;
    enm2 = na - 1;
   } //End if (incrFlag)
   else
    incrFlag = 1;
   
    /*Look for single small sub-diagonal element for l = en step -1 until low */

   for (var i = low; i <= en; i++){
    l = en + low - i;
    if (l == low)
     break;
    s = Math.abs(A[l - 1][l - 1]) + Math.abs(A[l][l]);
    if (s == 0.0)
     s = norm;
    tst1 = s;
    tst2 = tst1 + Math.abs(A[l][l - 1]);
    if (tst2 == tst1)
     break;
   } //End for i

   x = A[en][en];

   if (l == en){  //One root found
    wr[en] = A[en][en] = x + t;
    wi[en] = 0.0;
    en--;
    continue;
   } //End if (l == en)

   y = A[na][na];
   w = A[en][na]*A[na][en];
   
   if (l == na){  //Two roots found
    p = (-x + y)/2;
    q = p*p + w;
    zz = Math.sqrt(Math.abs(q));
    x = A[en][en] = x + t;
    A[na][na] = y + t;
    if (q >= 0.0){//Real Pair
     zz = ((p < 0.0) ? (-zz + p) : (p + zz));
     wr[en] = wr[na] = x + zz;
     if (zz != 0.0)
      wr[en] = -(w/zz) + x;
     wi[en] = wi[na] = 0.0;
     x = A[en][na];
     s = Math.abs(x) + Math.abs(zz);
     p = x/s;
     q = zz/s;
     r = Math.sqrt(p*p + q*q);
     p /= r;
     q /= r;
     for (var j = na; j < N; j++){ //Row modification
      zz = A[na][j];
      A[na][j] = q*zz + p*A[en][j];
      A[en][j] = -(p*zz) + q*A[en][j];
     }//End for j
     for (var j = 0; j <= en; j++){ // Column modification
      zz = A[j][na];
      A[j][na] = q*zz + p*A[j][en];
      A[j][en] = -(p*zz) + q*A[j][en];
     }//End for j
     for (var j = low; j <= igh; j++){//Accumulate transformations
      zz = B[j][na];
      B[j][na] = q*zz + p*B[j][en];
      B[j][en] = -(p*zz) + q*B[j][en];
     }//End for j
    } //End if (q >= 0.0)
    else {//else q < 0.0
     wr[en] = wr[na] = x + p;
     wi[na] = zz;
     wi[en] = -zz;
    } //End else
    en--;
    en--;
    continue;
   }//End if (l == na)
   
   if (itn == 0){ //Set error; all eigenvalues have not converged after 30 iterations.
    oPar.outEr = en + 1;
    return;
   }//End if (itn == 0)

   if ((its == 10) || (its == 20)){ //Form exceptional shift
    t += x;
    for (var i = low; i <= en; i++)
     A[i][i] += -x;
    s = Math.abs(A[en][na]) + Math.abs(A[na][enm2]);
    y = x = 0.75*s;
    w = -0.4375*s*s;
   } //End if (its equals 10 or 20)

   its++;
   itn--;

  /*Look for two consecutive small sub-diagonal elements. Do m = en - 2 to l in increments of -1 */

   for (m = enm2; m >= l; m--){
    zz = A[m][m];
    r = -zz + x;
    s = -zz + y;
    p = (-w + r*s)/A[m + 1][m] + A[m][m + 1];
    q = -(zz + r + s) + A[m + 1][m + 1] ;
    r = A[m + 2][m + 1];
    s = Math.abs(p) + Math.abs(q) + Math.abs(r);
    p /= s;
    q /= s;
    r /= s;
    if (m == l)
     break;
    tst1 = Math.abs(p) * (Math.abs(A[m - 1][m - 1]) + Math.abs(zz) + Math.abs(A[m + 1][m + 1]));
    tst2 = tst1 + Math.abs(A[m][m - 1]) * (Math.abs(q) + Math.abs(r));
    if (tst1 == tst2)
     break;
   }//End for m

   mp2 = m + 2;
   
   for (var i = mp2; i <= en; i++){
    A[i][i - 2] = 0.0;
    if (i == mp2)
     continue;
    A[i][i - 3] = 0.0;
   }//End for i
   
   /* Double qr step involving rows l to en and columns m to en. */

   for (var i = m; i <= na; i++){
    notlas = ((i != na) ? 1 : 0);
    if (i != m){
     p = A[i][i - 1];
     q = A[i + 1][i - 1];
     r = ((notlas) ? A[i + 2][i - 1] : 0.0);
     x = Math.abs(p) + Math.abs(q) + Math.abs(r);
     if (x == 0.0)      //Drop through rest of for i loop
      continue;
     p /= x;
     q /= x;
     r /= x;
    } //End if (i != m)

    s = Math.sqrt(p*p + q*q + r*r);
    if (p < 0.0)
     s = -s;

    if (i != m)
     A[i][i - 1] = -(s*x);
    else {
     if (l != m)
      A[i][i - 1] = -A[i][i - 1];
    }

    p += s;
    x = p/s;
    y = q/s;
    zz = r/s;
    q /= p;
    r /= p;
    k = ((i + 3 < en) ? i + 3 : en);

    if (notlas){ //Do row modification
     for (var j = i; j < N; j++) {
      p = A[i][j] + q*A[i + 1][j] + r*A[i + 2][j];
      A[i][j] += -(p*x);
      A[i + 1][j] += -(p*y);
      A[i + 2][j] += -(p*zz);
     }//End for j

     for (var j = 0; j <= k; j++) {//Do column modification
      p = x*A[j][i] + y*A[j][i + 1] + zz*A[j][i + 2];
      A[j][i] += -p;
      A[j][i + 1] += -(p*q);
      A[j][i + 2] += -(p*r);
     }//End for j
     
     for (var j = low; j <= igh; j++) {//Accumulate transformations
      p = x*B[j][i] + y*B[j][i + 1] + zz*B[j][i + 2];
      B[j][i] += -p;
      B[j][i + 1] += -(p*q);
      B[j][i + 2] += -(p*r);
     } // End for j
    }//End if notlas

    else {
     for (var j = i; j < N; j++) {//Row modification
      p = A[i][j] + q*A[i + 1][j];
      A[i][j] += -(p*x);
      A[i + 1][j] += -(p*y);
     }//End for j

     for (var j = 0; j <= k; j++){//Column modification
      p = x*A[j][i] + y*A[j][i +1];
      A[j][i] += -p;
      A[j][i + 1] += -(p*q);
     }//End for j

     for (var j = low; j <= igh; j++){//Accumulate transformations
      p = x*B[j][i] + y*B[j][i +1];
      B[j][i] += -p;
      B[j][i + 1] += -(p*q);
     }//End for j

    } //End else if notlas
   }//End for i
   incrFlag = 0;
  }//End while (en >= low)

  if (norm == 0.0)
   return;

  //Step from (N - 1) to 0 in steps of -1

  for (en = (N - 1); en >= 0; en--){
   p = wr[en];
   q = wi[en];
   na = en - 1;

   if (q > 0.0)
    continue;

   if (q == 0.0){//Real vector
    m = en;
    A[en][en] = 1.0;

    for (var j = na; j >= 0; j--){
     w = -p + A[j][j];
     r = 0.0;
     for (var ii = m; ii <= en; ii++)
      r += A[j][ii]*A[ii][en];

     if (wi[j] < 0.0){
      zz = w;
      s = r;
     }//End wi[j] < 0.0

     else {//wi[j] >= 0.0
      m = j;
      if (wi[j] == 0.0){
       t = w;
       if (t == 0.0){
        t = tst1 = norm;
        do {
         t *= 0.01;
         tst2 = norm + t;
        } while (tst2 > tst1);
       } //End if t == 0.0
       A[j][en] = -(r/t);
      }//End if wi[j] == 0.0

      else { //else wi[j] > 0.0; Solve real equations
       x = A[j][j + 1];
       y = A[j + 1][j];
       q = (-p + wr[j])*(-p + wr[j]) + wi[j]*wi[j];
       A[j][en] = t = (-(zz*r) + x*s)/q;
       A[j + 1][en] = ((Math.abs(x) > Math.abs(zz)) ? -(r + w*t)/x : -(s + y*t)/zz);
      }//End  else wi[j] > 0.0

      // Overflow control
      t = Math.abs(A[j][en]);
      if (t == 0.0)
       continue; //go up to top of for j loop
      tst1 = t;
      tst2 = tst1 + 1.0/tst1;
      if (tst2 > tst1)
       continue; //go up to top of for j loop
      for (var ii = j; ii <= en; ii++)
       A[ii][en] /= t;

     }//End else wi[j] >= 0.0

    }//End for j
    
   }      //End q == 0.0

   else {//else q < 0.0, complex vector
   //Last vector component chosen imaginary so that eigenvector matrix is triangular
   m = na;

   if (Math.abs(A[en][na]) <= Math.abs(A[na][en]))
    cdivA(0.0, -A[na][en], -p + A[na][na], q, A, na, na, en);
   else {
    A[na][na] = q/A[en][na];
    A[na][en] = -(-p + A[en][en])/A[en][na];
   } //End else (Math.abs(A[en][na] > Math.abs(A[na][en])

   A[en][na] = 0.0;
   A[en][en] = 1.0;

   for (var j = (na - 1); j >= 0; j--) {
    w = -p + A[j][j];
    sa = ra = 0.0;
    
    for (var ii = m; ii <= en; ii++) {
     ra += A[j][ii]*A[ii][na];
     sa += A[j][ii]*A[ii][en];
    } //End for ii

    if (wi[j] < 0.0){
     zz = w;
     r = ra;
     s = sa;
     continue;
    } //End if (wi[j] < 0.0)

    //else wi[j] >= 0.0
    m = j;
    if (wi[j] == 0.0)
     cdivA(-ra, -sa, w, q, A, j, na, en);
    else {// else wi[j] > 0.0; solve complex equations
     x = A[j][j + 1];
     y = A[j + 1][j];
     vr = -(q*q) + (-p + wr[j])*(-p + wr[j]) + wi[j]*wi[j];
     vi = (-p + wr[j])*2.0*q;
     if ((vr == 0.0) && (vi == 0.0)){
      tst1 = norm*(Math.abs(w) + Math.abs(q) + Math.abs(x) + Math.abs(y) + Math.abs(zz));
      vr = tst1;
      do {
       vr *= 0.01;
       tst2 = tst1 + vr;
      } while (tst2 > tst1);
     } //End if vr and vi == 0.0
     cdivA(-(zz*ra) + x*r + q*sa, -(zz*sa + q*ra) + x*s, vr, vi, A, j, na, en);

     if (Math.abs(x) > (Math.abs(zz) + Math.abs(q))){
      A[j + 1][na] = (-(ra + w*A[j][na]) + q*A[j][en])/x;
      A[j + 1][en] = -(sa + w*A[j][en] + q*A[j][na])/x;
     }//End if
     else
      cdivA(-(r + y*A[j][na]), -(s + y*A[j][en]), zz, q, A, j + 1, na, en);

    }//End else wi[j] > 0.0
      
    t = ((Math.abs(A[j][na]) >= Math.abs(A[j][en])) ? Math.abs(A[j][na]) : Math.abs(A[j][en]));

    if (t == 0.0)
     continue; // go to top of for j loop

    tst1 = t;
    tst2 = tst1 + 1.0/tst1;
    if (tst2 > tst1)
     continue; //go to top of for j loop
      
    for (var ii = j; ii <= en; ii++){
     A[ii][na] /= t;
     A[ii][en] /= t;
    } //End for ii loop

   } // End for j
    
   }//End else q < 0.0
   
  }//End for en

  //End back substitution. Vectors of isolated roots.

  for (var i = 0; i < N; i++){
   if ((i < low) || (i > igh)) {
    for (var j = i; j < N; j++)
     B[i][j] = A[i][j];
   }//End if i
  }//End for i

  // Multiply by transformation matrix to give vectors of original full matrix.

  //Step from (N - 1) to low in steps of -1.

  for (var i = (N - 1); i >= low; i--){

   m = ((i < igh) ? i : igh);
   
   for (var ii = low; ii <= igh; ii++){
    zz = 0.0;
    for (var jj = low; jj <= m; jj++)
     zz += B[ii][jj]*A[jj][i];
    B[ii][i] = zz;
   }//End for ii
  }//End of for i loop

  return;
  } //End of function hqr2

  function norVecC(N, Z, wi){// Normalizes the eigenvectors

  // This subroutine is based on the LINPACK routine SNRM2, written 25 October 1982, modified
  // on 14 October 1993 by Sven Hammarling of NAG Ltd.
  // I have further modified it for use in this Javascript routine, for use with a column
  // of an array rather than a column vector.
  //
  // Z - eigenvector Matrix
  // wi - eigenvalue vector

   var scale, ssq, absxi, dummy, norm;
    
   for (var j = 0; j < N; j++){ //Go through the columns of the vector array 
    scale = 0.0;
    ssq = 1.0;

    for (var i = 0; i < N; i++){
     if (Z[i][j] != 0){
      absxi = Math.abs(Z[i][j]);
      dummy = scale/absxi;
      if (scale < absxi){
       ssq = 1.0 + ssq*dummy*dummy;
       scale = absxi;
      }//End if (scale < absxi)
      else
       ssq += 1.0/dummy/dummy;
     }//End if (Z[i][j] != 0)
    } //End for i

    if (wi[j] != 0){// If complex eigenvalue, take into account imaginary part of eigenvector
     for (var i = 0; i < N; i++){
      if (Z[i][j + 1] != 0){
       absxi = Math.abs(Z[i][j + 1]);
       dummy = scale/absxi;
       if (scale < absxi){
        ssq = 1.0 + ssq*dummy*dummy;
        scale = absxi;
       }//End if (scale < absxi)
       else
        ssq += 1.0/dummy/dummy;
       }//End if (Z[i][j + 1] != 0)
     } //End for i
    }//End if (wi[j] != 0)

    norm = scale*Math.sqrt(ssq); //This is the norm of the (possibly complex) vector

    for (var i = 0; i < N; i++)
     Z[i][j] /= norm;

    if (wi[j] != 0){// If complex eigenvalue, also scale imaginary part of eigenvector
     j++;
     for (var i = 0; i < N; i++)
     Z[i][j] /= norm;
    }//End if (wi[j] != 0)

   }// End for j
   
   return;
  } // End norVecC

  function calcEigSysReal(N, A, B, wr, wi, oPar){

  var scale = new Array(N);    //Contains information about transformations.
  var trace = new Array(N);    //Records row and column interchanges
  var radix = 2;               //Base of machine floating point representation.
  var c, f, g, r, s, b2 = radix*radix;
  var ierr = -1, igh, low, k = 0, l = N - 1, noconv;

  /* Balance the matrix to improve accuracy of eigenvalues. Introduces no rounding errors, since it scales A by powers of the radix.
  */

  //Search through rows, isolating eigenvalues and pushing them down.

  noconv = l;

  while (noconv >= 0){
   r = 0;

   for (var j = 0; j <= l; j++) {
    if (j == noconv) continue;
    if (A[noconv][j] != 0.0){
     r = 1;
     break;
    }
   } //End for j

   if (r == 0){
    scale[l] = noconv;

    if (noconv != l){
     for (var i = 0; i <= l; i++){
      f = A[i][noconv];
      A[i][noconv] = A[i][l];
      A[i][l] = f;
     }//End for i
     for (var i = 0; i < N; i++){
      f = A[noconv][i];
      A[noconv][i] = A[l][i];
      A[l][i] = f;
     }//End for i
    }//End if (noconv != l)

    if (l == 0)
     break;  //break out of while loop

    noconv = --l;

   }//End if (r == 0)

   else //else (r != 0)
    noconv--;

  }//End while (noconv >= 0)

  if (l > 0) {  //Search through columns, isolating eigenvalues and pushing them left.

   noconv = 0;

   while (noconv <= l){  
    c = 0;

    for (var i = k; i <= l; i++){
     if (i == noconv) continue;
     if (A[i][noconv] != 0.0){
      c = 1;
      break;
     }
    }//End for i

    if (c == 0){
     scale[k] = noconv;

     if (noconv != k){
      for (var i = 0; i <= l; i++){
       f = A[i][noconv];
       A[i][noconv] = A[i][k];
       A[i][k] = f;
      }//End for i
      for (var i = k; i < N; i++){
       f = A[noconv][i];
       A[noconv][i] = A[k][i];
       A[k][i] = f;
      }//End for i

     }//End if (noconv != k)

     noconv = ++k;
    }//End if (c == 0)

    else  //else (c != 0)
     noconv++;

   }//End while noconv

   //Balance the submatrix in rows k through l.

   for (var i = k; i <= l; i++)
    scale[i] = 1.0;

   //Iterative loop for norm reduction
   do {
    noconv = 0;
    for (var i = k; i <= l; i++) {
     c = r = 0.0;
     for (var j = k; j <= l; j++){
      if (j == i) continue;
      c += Math.abs(A[j][i]);
      r += Math.abs(A[i][j]);
     } // End for j
     if ((c == 0.0) || (r == 0.0)) continue;   //guard against zero c or r due to underflow
     g = r/radix;
     f = 1.0;
     s = c + r;
     while (c < g) {
      f *= radix;
      c *= b2;
     } // End while (c < g)
     g = r*radix;
     while (c >= g) {
      f /= radix;
      c /= b2;
     } // End while (c >= g)

     //Now balance
     if ((c + r)/f < 0.95*s) {
      g = 1.0/f;
      scale[i] *= f;
      noconv = 1;
      for (var j = k; j < N; j++)
       A[i][j] *= g;
      for (var j = 0; j <= l; j++)
       A[j][i] *= f;
     } //End if ((c + r)/f < 0.95*s)
    } // End for i
   } while (noconv);  // End of do-while loop.

  } //End if (l > 0)

  low = k;
  igh = l;

  //End of balanc
   
  /* Now reduce the real general Matrix to upper-Hessenberg form using stabilized elementary similarity transformations. */

  for (var i = (low + 1); i < igh; i++){
   k = i;
   c = 0.0;

   for (var j = i; j <= igh; j++){
    if (Math.abs(A[j][i - 1]) > Math.abs(c)){
     c = A[j][i - 1];
     k = j;
    }//End if
   }//End for j

   trace[i] = k;

   if (k != i){
    for (var j = (i - 1); j < N; j++){
     r = A[k][j];
     A[k][j] = A[i][j];
     A[i][j] = r;
    }//End for j

    for (var j = 0; j <= igh; j++){
     r = A[j][k];
     A[j][k] = A[j][i];
     A[j][i] = r;
    }//End for j
   }//End if (k != i)

   if (c != 0.0){
    for (var m = (i + 1); m <= igh; m++){
     r = A[m][i - 1];

     if (r != 0.0){
      r = A[m][i - 1] = r/c;
      for (var j = i; j < N; j++)
       A[m][j] += -(r*A[i][j]);
      for (var j = 0; j <= igh; j++)
       A[j][i] += r*A[j][m];
     }//End if (r != 0.0)
    }//End for m
   }//End if (c != 0)
  }  //End for i.

  /* Accumulate the stabilized elementary similarity transformations used in the reduction of A to upper-Hessenberg form. Introduces no rounding errors since it only transfers the multipliers used in the reduction process into the eigenvector matrix. */

  for (var i = 0; i < N; i++){ //Initialize B to the identity Matrix.
   for (var j = 0; j < N; j++)
    B[i][j] = 0.0;
   B[i][i] = 1.0;
  } //End for i

  for (var i = (igh - 1); i >= (low + 1); i--){
   k = trace[i];
   for (var j = (i + 1); j <= igh; j++)
    B[j][i] = A[j][i - 1];

   if (i == k)
    continue;

   for (var j = i; j <= igh; j++){
    B[i][j] = B[k][j];
    B[k][j] = 0.0;
   } //End for j

   B[k][i] = 1.0; 

  } // End for i
   
  oPar.outEr = ierr;
  hqr2(N, A, B, low, igh, wi, wr, oPar);
  ierr = oPar.outEr;

  if (ierr == -1){

   if (low != igh){
    for (var i = low; i <= igh; i++){
     s = scale[i];
     for (var j = 0; j < N; j++)
      B[i][j] *= s;
    }//End for i
   }//End if (low != igh)

   for (var i = (low - 1); i >= 0; i--){
    k = scale[i];
    if (k != i){
     for (var j = 0; j < N; j++){
      s = B[i][j];
      B[i][j] = B[k][j];
      B[k][j] = s;
     }//End for j
    }//End if k != i
   }//End for i

   for (var i = (igh + 1); i < N; i++){
    k = scale[i];
    if (k != i){
     for (var j = 0; j < N; j++){
      s = B[i][j];
      B[i][j] = B[k][j];
      B[k][j] = s;
     }//End for j
    }//End if k != i
   }//End for i

   norVecC(N, B, wi);  //Normalize the eigenvectors
    
   }//End if ierr = -1

   return;
  }  //End of compRealEigSys


  exports.calcEigr = function(matrix){
    // Main program for computing the eigenvalues and eigenvectors of a square, real matrix, A.
    // Assumes all entries are real numbers.

   var MATRIXDIM = matrix.length;	// The dimension, N, of the system to be entered.

   var B = new Array(MATRIXDIM);

   var wi = new Array(MATRIXDIM);
   var wr = new Array(MATRIXDIM);

   for (var i = 0; i < MATRIXDIM; i++){	// Create the square Matrices
     B[i] = new Array(MATRIXDIM);
     } // End of for i loop

   // outputPar is a dummy variable for passing output parameters by reference
   var outputPar = new Object();
   outputPar.outEr = 0;			// Error code, ierr

   // document.getElementById('output').innerHTML = "";

   //  *******************************************************************
   // At this point, MATRIXDIM should be the matrix dimension
   // and the A matrix should contain its appropriate values.
   //  ********************************************************************
   calcEigSysReal(MATRIXDIM, matrix, B, wr, wi, outputPar);


  return {matrix: B, wi: wi, wr: wr};
  }

})();
