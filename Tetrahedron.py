import sympy as sp

################################################################################
#
def analyze_tetrahedron_exact(AB, AC, AD, BC, BD, CD):
  """Computes exact symbolic geometric properties, canonical 3D vertex coordinates,

  and special centers (centroid, incenter, circumcenter) of a tetrahedron.

  Inputs can be integers, floats (automatically converted to exact rationals),
  or SymPy expressions (e.g., sp.sqrt(2)).
  """
  
  # 0. Exact Edge Lengths ######################################################

  AB, AC, AD, BC, BD, CD = map(sp.nsimplify, (AB, AC, AD, BC, BD, CD))

  edge = {"Edge AB": sp.simplify(AB),
          "Edge AC": sp.simplify(AC),
          "Edge AD": sp.simplify(AD),
          "Edge BC": sp.simplify(BC),
          "Edge BD": sp.simplify(BD),
          "Edge CD": sp.simplify(CD)}

  # 1. Exact Face Areas ########################################################

  def exact_triangle_area(x, y, z):
    s = sp.simplify((x + y + z) / 2)
    return sp.simplify(sp.sqrt(s * (s - x) * (s - y) * (s - z)))

  ABC = exact_triangle_area(AB, AC, BC)
  ABD = exact_triangle_area(AB, AD, BD)
  ACD = exact_triangle_area(AC, AD, CD)
  BCD = exact_triangle_area(BC, BD, CD)

  area = {"Area ABC": ABC,
          "Area ABD": ABD,
          "Area ACD": ACD,
          "Area BCD": BCD}
  
  total_area = sp.simplify(ABC + ABD + ACD + BCD)  

  # 2. Exact Volume via Piero della Francesca's / Cayley-Menger formula ########

  AB2, AC2, AD2 = map(sp.nsimplify, (AB**2, AC**2, AD**2))
  BC2, BD2, CD2 = map(sp.nsimplify, (BC**2, BD**2, CD**2))
  
  volume = sp.simplify(sp.sqrt(BC2 * AD2 * (AC2 + AB2 + BD2 + CD2 - BC2 - AD2)
                             + AC2 * BD2 * (BC2 + AB2 + AD2 + CD2 - AC2 - BD2)
                             + AB2 * CD2 * (BC2 + AC2 + AD2 + BD2 - AB2 - CD2)
                             - BC2 * AC2 * AB2 - BC2 * BD2 * CD2
                             - AC2 * AD2 * CD2 - AB2 * AD2 * BD2)/12)

  # 3. Canonical 3D Vertex Coordinates #########################################
  
  Cx_Dx = sp.simplify((AB**2 + AC**2 - BC**2)*(AB**2 + AD**2 - BD**2)/(2*AB)**2) 
  Cy    = sp.simplify((2*ABC) / AB)
  
  A = (sp.Integer(0), sp.Integer(0), sp.Integer(0))
  B = (AB,            sp.Integer(0), sp.Integer(0))
  C = (sp.simplify((AB**2 + AC**2 - BC**2) / (2*AB)), 
       sp.simplify((2*ABC) / AB), sp.Integer(0))
  D = (sp.simplify((AB**2 + AD**2 - BD**2) / (2*AB)),
       sp.simplify((AC**2 + AD**2 - CD**2 - 2* Cx_Dx) / (2 * Cy)),
       sp.simplify((3 * volume) / ABC))
  
  vertex = {"Vertex A": A, 
            "Vertex B": B, 
            "Vertex C": C, 
            "Vertex D": D}

  # 4. Special Centers in 3D Space #############################################

  Oy = sp.simplify((C[0]**2 + C[1]**2 - AB * C[0]) / (2*C[1]))
  centroid = tuple(sp.simplify((A[i] + B[i] + C[i] + D[i])/4) for i in range(3))
  incenter = tuple(sp.simplify((BCD*A[i] + ACD*B[i] + ABD*C[i] + ABC*D[i]) 
                              / total_area) for i in range(3))
  circumcenter = (sp.simplify(AB/2),
                  sp.simplify((C[0]**2 + C[1]**2 - AB * C[0]) / (2*C[1])), 
                  sp.simplify((D[0]**2 + D[1]**2 + D[2]**2 
                             - AB*D[0] - 2*D[1]*Oy) / (2*D[2])))

  centers = {"Centroid":     centroid,
             "Incenter":     incenter,
             "Circumcenter": circumcenter}

  # 5. Exact Inradius and Circumradius #########################################
  
  inradius     = sp.simplify((3*volume) / total_area)
  circumradius = sp.simplify(exact_triangle_area(AB*CD,AC*BD,AD*BC)/(6*volume))

  # 6. Exact Face Heights / Altitudes ##########################################
  
  heights = {"hA (from A to BCD)": sp.simplify((3 * volume) / BCD),
             "hB (from B to ACD)": sp.simplify((3 * volume) / ACD),
             "hC (from C to ABD)": sp.simplify((3 * volume) / ABD),
             "hD (from D to ABC)": sp.simplify((3 * volume) / ABC)}

  # 7. Exact 12 Face Angles (in radians) #######################################

  angles = {"Face ABC": {
              "angle A": sp.simplify(sp.acos((AB**2+AC**2-BC**2)/(2*AB*AC))),
              "angle B": sp.simplify(sp.acos((AB**2+BC**2-AC**2)/(2*AB*BC))),
              "angle C": sp.simplify(sp.acos((AC**2+BC**2-AB**2)/(2*AC*BC)))},
            "Face ABD": {
              "angle A": sp.simplify(sp.acos((AB**2+AD**2-BD**2)/(2*AB*AD))),
              "angle B": sp.simplify(sp.acos((AB**2+BD**2-AD**2)/(2*AB*BD))),
              "angle D": sp.simplify(sp.acos((AD**2+BD**2-AB**2)/(2*AD*BD)))},
            "Face ACD": {
              "angle A": sp.simplify(sp.acos((AC**2+AD**2-CD**2)/(2*AC*AD))),
              "angle C": sp.simplify(sp.acos((AC**2+CD**2-AD**2)/(2*AC*CD))),
              "angle D": sp.simplify(sp.acos((AD**2+CD**2-AC**2)/(2*AD*CD)))},
            "Face BCD": {
              "angle B": sp.simplify(sp.acos((BC**2+BD**2-CD**2)/(2*BC*BD))),
              "angle C": sp.simplify(sp.acos((BC**2+CD**2-BD**2)/(2*BC*CD))),
              "angle D": sp.simplify(sp.acos((BD**2+CD**2-BC**2)/(2*BD*CD)))}}

  # 8. Exact 6 Dihedral Angles (in radians) ####################################
  
  def exact_dihedral(AB, AC, AD, BC, BD, CD):
    """
    Returns dihedral angle of edge AB using the spherical law of cosines at A.
    """
    
    cos_BAC = sp.simplify((AB**2 + AC**2 - BC**2) / (2*AB*AC))
    cos_BAD = sp.simplify((AB**2 + AD**2 - BD**2) / (2*AB*AD))
    cos_CAD = sp.simplify((AC**2 + AD**2 - CD**2) / (2*AC*AD))
    
    sin_BAC = sp.simplify(sp.sqrt(1 - cos_BAC**2))
    sin_BAD = sp.simplify(sp.sqrt(1 - cos_BAD**2))
    
    num = sp.simplify(cos_CAD - cos_BAC * cos_BAD)
    den = sp.simplify(sin_BAC*sin_BAD)
    
    return sp.simplify(sp.acos(num/den))
    
  dihedral = {"AB": exact_dihedral(AB, AC, AD,   BC, BD, CD),
              "AC": exact_dihedral(AC, AD, AB,   CD, BC, BD),
              "AD": exact_dihedral(AD, AB, AC,   BD, CD, BC),
              "BC": exact_dihedral(BC, AB, BD,   AC, CD, AD),
              "BD": exact_dihedral(BD, BC, AB,   CD, AD, AC),
              "CD": exact_dihedral(CD, AD, BD,   AC, BC, AB)}

  ##############################################################################

  return {"vertex":       vertex,
          "edge":         edge,
          "centers":      centers,
          "volume":       volume,
          "total_area":   total_area,
          "inradius":     inradius,
          "circumradius": circumradius,
          "area":         area,
          "heights":      heights,
          "angles":       angles,
          "dihedral":     dihedral}
#
################################################################################
#
def print_exact_summary(props):
  """Prints a clean summary showing exact symbolic values alongside decimals."""

  def fmt_val(expr):
    return f"{expr}   (≈ {float(expr.evalf()):.4f})"

  def fmt_ang(expr):
    deg = float((expr * 180 / sp.pi).evalf())
    return f"{expr} rad   (≈ {deg:.2f}°)"

  def fmt_pt(pt):
    exact_str = f"({pt[0]}, {pt[1]}, {pt[2]})"
    approx_str = (f"({float(pt[0].evalf()):.4f},"
                  f" {float(pt[1].evalf()):.4f},"
                  f" {float(pt[2].evalf()):.4f})")
    return f"{exact_str}  (≈ {approx_str})"

  ##############################################################################

  print("==========================================================")
  print("                 EXACT TETRAHEDRON SUMMARY                ")
  print("==========================================================")
  
  print("\n--- GENERAL METRIC PROPERTIES ---\n")
  print(f"  Volume       : {fmt_val(props['volume'])}")
  print(f"  Surface Area : {fmt_val(props['total_area'])}")
  print(f"  Inradius     : {fmt_val(props['inradius'])}")
  print(f"  Circumradius : {fmt_val(props['circumradius'])}")

  print("\n--- CANONICAL 3D VERTEX COORDINATES ---\n")
  for v, pt in props["vertex"].items():
    print(f"  {v:2}: {fmt_pt(pt)}")

  print("\n--- 3D SPECIAL CENTERS ---\n")
  for name, pt in props["centers"].items():
    print(f"  {name:13}: {fmt_pt(pt)}")

  print("\n--- FACE AREAS ---\n")
  for face, area in props["area"].items():
    print(f"  {face:17}: {fmt_val(area)}")

  print("\n--- FACE HEIGHTS (ALTITUDES) ---\n")
  for alt, h in props["heights"].items():
    print(f"  {alt:19}: {fmt_val(h)}")

  print("\n--- FACE ANGLES ---\n")
  for face, angles in props["angles"].items():
    print(f"  {face}:")
    for name, ang in angles.items():
      print(f"    {name:14}: {fmt_ang(ang)}")

  print("\n--- DIHEDRAL ANGLES ---\n")
  for edge, ang in props["dihedral"].items():
    print(f"  {edge:3}: {fmt_ang(ang)}")

  print("\n--- EDGE LENGTHS ---\n")
  for edge, length in props["edge"].items():
    print(f"  {edge:3}: {fmt_val(length)}")
#
################################################################################
#
if __name__ == "__main__":
  
  from math import sqrt

  COS_18 = sqrt(5/8 + sqrt(5)/8)
  COS_30 = sqrt(3)/2
  PHI    = (1+sqrt(5))/2
  
  B0 = 1
  B1 = PHI
  B2 = (PHI+1)
  
  Y0 = COS_30*1
  Y1 = COS_30*PHI
  Y2 = COS_30*(PHI+1)
  
  R0 = COS_18*1
  R1 = COS_18*PHI
  R2 = COS_18*(PHI+1)
  
  G0 = sqrt(2)*1
  G1 = sqrt(2)*PHI
  G2 = sqrt(2)*(PHI+1)
  
  ZOME = {56:(B0,Y0,Y0,Y0,Y0,B0)}
  
  props = analyze_tetrahedron_exact(*ZOME[56])
  
  #props = analyze_tetrahedron_exact(1,1,1,1,1,1)  
  #props = analyze_tetrahedron_exact(2,sqrt(3),sqrt(3),sqrt(3),sqrt(3),2)
  #props = analyze_tetrahedron_exact(2,2,2,sqrt(3),sqrt(3),sqrt(3))
  #props = analyze_tetrahedron_exact(1,1,2,sqrt(2),sqrt(5),sqrt(5))
  
  print_exact_summary(props)
#
################################################################################
