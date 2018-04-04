/*
 * Implement an abstract data type "Quaternion" according to http://enos.itcollege.ee/%7Ejpoial/algoritmid/quat/
 * Used materials:
 * * http://opik.fyysika.ee/index.php/book/section/4592
 * * https://en.wikipedia.org/wiki/Quaternion
 * * https://www.youtube.com/watch?v=jlskQDR8-bY
 * * http://enos.itcollege.ee/~jpoial/algorithms/examples/Num.java
 * * https://stackoverflow.com/questions/2808535/round-a-double-to-2-decimal-places > https://stackoverflow.com/a/25645952
 * * https://stackoverflow.com/questions/40618710/regex-to-exclude-first-and-last-characters-of-match
 * * https://stackoverflow.com/questions/364454/findbugs-warning-equals-method-should-not-assume-anything-about-the-type-of-its
 *
 * Kommentaar: Ümardamise teema ei tohiks olla tehetes oluline - ainult võrdluses ja pisut ka sisendis/väljundis. Samas
 * isZero peaks kasutama ümardamist, meetodid inverse ja divideBy... peaksid kontrollima nulliga jagamist.
 */

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;

/** Quaternions. Basic operations. */
public class Quaternion {
   /** Constructor from four double values.
    * @param a real part
    * @param b imaginary part i
    * @param c imaginary part j
    * @param d imaginary part k
    */
   private double rPart;
   private double iPart;
   private double jPart;
   private double kPart;
   public Quaternion (double a, double b, double c, double d) {
      this.rPart = a;
      this.iPart = b;
      this.jPart = c;
      this.kPart = d;
   }

   /** Rounded double.
    * @return rounded double
    */
   private double roundDouble(double d){
      return new BigDecimal(((Double)d).toString()).setScale(6, RoundingMode.HALF_UP).doubleValue();
   }

   /** Real part of the quaternion.
    * @return real part
    */
   public double getRpart() {
      return this.rPart;
   }

   /** Imaginary part i of the quaternion. 
    * @return imaginary part i
    */
   public double getIpart() {
      return this.iPart;
   }

   /** Imaginary part j of the quaternion. 
    * @return imaginary part j
    */
   public double getJpart() {
      return this.jPart;
   }

   /** Imaginary part k of the quaternion. 
    * @return imaginary part k
    */
   public double getKpart() {
      return this.kPart;
   }

   /** Conversion of the quaternion to the string.
    * @return a string form of this quaternion: 
    * "a+bi+cj+dk"
    * (without any brackets)
    */
   @Override
   public String toString() {
      String s = String.valueOf(rPart) + "+" + String.valueOf(iPart) + "i+" + String.valueOf(jPart) + "j+"
              + String.valueOf(kPart) + "k";

      return s.replace("+-", "-");
   }

   /** Conversion from the string to the quaternion. 
    * Reverse to <code>toString</code> method.
    * @throws IllegalArgumentException if string s does not represent 
    *     a quaternion (defined by the <code>toString</code> method)
    * @param s string of form produced by the <code>toString</code> method
    * @return a quaternion represented by string s
    */
   public static Quaternion valueOf (String s) {
      if (s.isEmpty()) throw new RuntimeException("Entered string \"" + s + "\" does not contain any information!");
      if (!s.contains("i")) throw new RuntimeException("Entered string \"" + s + "\" does not contain imaginary part \"i\"!");
      if (!s.contains("j")) throw new RuntimeException("Entered string \"" + s + "\" does not contain imaginary part \"j\"!");
      if (!s.contains("k")) throw new RuntimeException("Entered string \"" + s + "\" does not contain imaginary part \"k\"!");
      String signReplacement = s.replaceAll("(?<=.)[-]+", "+-");
      String[] stringSplit = signReplacement.split("[+]++");
      if (stringSplit.length < 4) throw new RuntimeException("Entered string \"" + s + "\" does not contain real part value!");
      if (stringSplit.length > 4) throw new RuntimeException("Entered string \"" + s + "\" contains too many values!");
      Double[] quaternionValues = new Double[4];
      for (String item: stringSplit) {
         int index = 0;
         if (item.contains("i")) index = 1;
         if (item.contains("j")) index = 2;
         if (item.contains("k")) index = 3;
         String cleanItem = item.replaceAll("[ijk]","");
         if (cleanItem.isEmpty()) cleanItem = "1";
         if (cleanItem.equals("-")) cleanItem = "-1";
         double value = Double.parseDouble(cleanItem);
         quaternionValues[index] = value;
      }
      double r = quaternionValues[0];
      double i = quaternionValues[1];
      double j = quaternionValues[2];
      double k = quaternionValues[3];
      return new Quaternion(r, i, j, k);
   }

   /** Clone of the quaternion.
    * @return independent clone of <code>this</code>
    */
   @Override
   public Object clone() throws CloneNotSupportedException {
      String string = this.toString();
      return Quaternion.valueOf(string);
   }

   /** Test whether the quaternion is zero. 
    * @return true, if the real part and all the imaginary parts are (close to) zero
    */
   public boolean isZero() {
      if (roundDouble(this.rPart) == 0) {
         if (roundDouble(this.iPart) == 0) {
            if (roundDouble(this.jPart) == 0) {
               if (roundDouble(this.kPart) == 0) {
                  return true;
               }
            }
         }
      }
      return false;
   }

   /** Conjugate of the quaternion. Expressed by the formula 
    *     conjugate(a+bi+cj+dk) = a-bi-cj-dk
    * @return conjugate of <code>this</code>
    */
   public Quaternion conjugate() {
      double conjIPart = -this.iPart;
      double conjJPart = -this.jPart;
      double conjKPart = -this.kPart;
      return new Quaternion(this.rPart, conjIPart, conjJPart, conjKPart);
   }

   /** Opposite of the quaternion. Expressed by the formula 
    *    opposite(a+bi+cj+dk) = -a-bi-cj-dk
    * @return quaternion <code>-this</code>
    */
   public Quaternion opposite() {
      double oppRPart = -this.rPart;
      double oppIPart = -this.iPart;
      double oppJPart = -this.jPart;
      double oppKPart = -this.kPart;
      return new Quaternion(oppRPart, oppIPart, oppJPart, oppKPart);
   }

   /** Sum of quaternions. Expressed by the formula 
    *    (a1+b1i+c1j+d1k) + (a2+b2i+c2j+d2k) = (a1+a2) + (b1+b2)i + (c1+c2)j + (d1+d2)k
    * @param q addend
    * @return quaternion <code>this+q</code>
    */
   public Quaternion plus (Quaternion q) {
      double plusRPart = this.rPart + q.rPart;
      double plusIPart = this.iPart + q.iPart;
      double plusJPart = this.jPart + q.jPart;
      double plusKPart = this.kPart + q.kPart;
      return new Quaternion(plusRPart, plusIPart, plusJPart, plusKPart);
   }

   /** Product of quaternions. Expressed by the formula
    *  (a1+b1i+c1j+d1k) * (a2+b2i+c2j+d2k) = (a1a2-b1b2-c1c2-d1d2) + (a1b2+b1a2+c1d2-d1c2)i +
    *  (a1c2-b1d2+c1a2+d1b2)j + (a1d2+b1c2-c1b2+d1a2)k
    * @param q factor
    * @return quaternion <code>this*q</code>
    */
   public Quaternion times (Quaternion q) {
      double a1 = this.rPart;
      double a2 = q.rPart;
      double b1 = this.iPart;
      double b2 = q.iPart;
      double c1 = this.jPart;
      double c2 = q.jPart;
      double d1 = this.kPart;
      double d2 = q.kPart;
      double qTimesRPart = roundDouble(a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2);
      double qTimesIPart = roundDouble(a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2);
      double qTimesJPart = roundDouble(a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2);
      double qTimesKPart = roundDouble(a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2);
      return new Quaternion(qTimesRPart, qTimesIPart, qTimesJPart, qTimesKPart);
   }

   /** Multiplication by a coefficient.
    * @param r coefficient
    * @return quaternion <code>this*r</code>
    */
   public Quaternion times (double r) {
      double rTimesRPart = roundDouble(r * this.rPart);
      double rTimesIPart = roundDouble(r * this.iPart);
      double rTimesJPart = roundDouble(r * this.jPart);
      double rTimesKPart = roundDouble(r * this.kPart);
      return new Quaternion(rTimesRPart, rTimesIPart, rTimesJPart, rTimesKPart);
   }

   /** Inverse of the quaternion. Expressed by the formula
    *     1/(a+bi+cj+dk) = a/(a*a+b*b+c*c+d*d) + 
    *     ((-b)/(a*a+b*b+c*c+d*d))i + ((-c)/(a*a+b*b+c*c+d*d))j + ((-d)/(a*a+b*b+c*c+d*d))k
    * @return quaternion <code>1/this</code>
    */
   public Quaternion inverse() {
      if (this.isZero()) throw new RuntimeException("Inverse is not possible! Quaternion \"" + this.toString() +"\"," +
              " is equal to zero!");
      double a = this.rPart;
      double b = this.iPart;
      double c = this.jPart;
      double d = this.kPart;
      double divideBy = a * a + b * b + c * c + d * d;
      double inverseA = roundDouble(a / divideBy);
      double inverseB = roundDouble(-b / divideBy);
      double inverseC = roundDouble(-c / divideBy);
      double inverseD = roundDouble(-d / divideBy);
      return new Quaternion(inverseA, inverseB, inverseC, inverseD);
   }

   /** Difference of quaternions. Expressed as addition to the opposite.
    * @param q subtrahend
    * @return quaternion <code>this-q</code>
    */
   public Quaternion minus (Quaternion q) {
      Quaternion oppositeQ = q.opposite();
      return this.plus(oppositeQ);
   }

   /** Right quotient of quaternions. Expressed as multiplication to the inverse.
    * @param q (right) divisor
    * @return quaternion <code>this*inverse(q)</code>
    */
   public Quaternion divideByRight (Quaternion q) {
      if (q.isZero()) throw new RuntimeException("Can not divide by given quaternion \"" + q.toString() + "\", it equals" +
              " to zero!");
      Quaternion inverseQ = q.inverse();
      return this.times(inverseQ);
   }

   /** Left quotient of quaternions.
    * @param q (left) divisor
    * @return quaternion <code>inverse(q)*this</code>
    */
   public Quaternion divideByLeft (Quaternion q) {
      if (q.isZero()) throw new RuntimeException("Can not divide the given quaternion \"" + q.toString() + "\". The" +
              " quaternion equals to zero!");
      Quaternion inverseQ = q.inverse();
      return inverseQ.times(this);
   }
   
   /** Equality test of quaternions. Difference of equal numbers
    *     is (close to) zero.
    * @param qo second quaternion
    * @return logical value of the expression <code>this.equals(qo)</code>
    */
   @Override
   public boolean equals (Object qo) {
      if (!(qo instanceof Quaternion)) return false;
      double qoRPart = ((Quaternion)qo).rPart;
      double qoIPart = ((Quaternion)qo).iPart;
      double qoJPart = ((Quaternion)qo).jPart;
      double qoKPart = ((Quaternion)qo).kPart;
      return (roundDouble(qoRPart) == roundDouble(this.rPart) && roundDouble(qoIPart) == roundDouble(this.iPart) &&
              roundDouble(qoJPart) == roundDouble(this.jPart) && roundDouble(qoKPart) == roundDouble(this.kPart));
   }

   /** Dot product of quaternions. (p*conjugate(q) + q*conjugate(p))/2
    * @param q factor
    * @return dot product of this and q
    */
   public Quaternion dotMult (Quaternion q) {
      Quaternion conjugateP = this.conjugate();
      Quaternion conjugateQ = q.conjugate();
      Quaternion pTimesConjQ = this.times(conjugateQ);
      Quaternion qTimesConjP = q.times(conjugateP);
      Quaternion pQPlusQP = pTimesConjQ.plus(qTimesConjP);
      return pQPlusQP.times(0.5);
   }

   /** Integer hashCode has to be the same for equal objects.
    * @return hashcode
    */
   @Override
   public int hashCode() {
      return this.toString().hashCode();
   }

   /** Norm of the quaternion. Expressed by the formula 
    *     norm(a+bi+cj+dk) = Math.sqrt(a*a+b*b+c*c+d*d)
    * @return norm of <code>this</code> (norm is a real number)
    */
   public double norm() {
      double sqrRPart = this.rPart * this.rPart;
      double sqrIPart = this.iPart * this.iPart;
      double sqrJPart = this.jPart * this.jPart;
      double sqrKPart = this.kPart * this.kPart;
      return Math.sqrt(sqrRPart + sqrIPart + sqrJPart + sqrKPart);
   }

   /** Main method for testing purposes. 
    * @param arg command line parameters
    */
   public static void main (String[] arg) {
      Quaternion arv1 = new Quaternion (-1., 1, 2., -2.);
      if (arg.length > 0)
         arv1 = valueOf (arg[0]);
      System.out.println ("first: " + arv1.toString());
      System.out.println ("real: " + arv1.getRpart());
      System.out.println ("imagi: " + arv1.getIpart());
      System.out.println ("imagj: " + arv1.getJpart());
      System.out.println ("imagk: " + arv1.getKpart());
      System.out.println ("isZero: " + arv1.isZero());
      System.out.println ("conjugate: " + arv1.conjugate());
      System.out.println ("opposite: " + arv1.opposite());
      System.out.println ("hashCode: " + arv1.hashCode());
      Quaternion res = null;
      try {
         res = (Quaternion)arv1.clone();
      } catch (CloneNotSupportedException e) {}
      System.out.println ("clone equals to original: " + res.equals (arv1));
      System.out.println ("clone is not the same object: " + (res!=arv1));
      System.out.println ("hashCode: " + res.hashCode());
      res = valueOf (arv1.toString());
      System.out.println ("string conversion equals to original: " 
         + res.equals (arv1));
      Quaternion arv2 = new Quaternion (1., -2.,  -1., 2.);
      if (arg.length > 1)
         arv2 = valueOf (arg[1]);
      System.out.println ("second: " + arv2.toString());
      System.out.println ("hashCode: " + arv2.hashCode());
      System.out.println ("equals: " + arv1.equals (arv2));
      res = arv1.plus (arv2);
      System.out.println ("plus: " + res);
      System.out.println ("times: " + arv1.times (arv2));
      System.out.println ("minus: " + arv1.minus (arv2));
      double mm = arv1.norm();
      System.out.println ("norm: " + mm);
      System.out.println ("inverse: " + arv1.inverse());
      System.out.println ("divideByRight: " + arv1.divideByRight (arv2));
      System.out.println ("divideByLeft: " + arv1.divideByLeft (arv2));
      System.out.println ("dotMult: " + arv1.dotMult (arv2));
      Quaternion arv3 = valueOf("0i+0j-0k-0");
      System.out.println(arv3.toString());
      arv3.isZero();
      Quaternion arv4 = arv3.divideByLeft(arv1);
      System.out.println(arv4);
      arv3.inverse();
   }
}
// end of file
