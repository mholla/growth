import sys

##### PRELIM ##########################################################

name = sys.argv[1]
orientation = sys.argv[2]
umat_filename = "../umat_" + name + ".f"
perturb = ['','','']

##### Cauchy_11-driven volumetric growth ##############################
# sdvini = "sdvini.f"
# umat = "callumat_cauchy.f"
# umat_cort = "../../6 Growth Models/Morphogenetic/Fiber/umat_length_linear.f"
# cor_start = 36
# umat_intro = "brain-cauchy.f"
# umat_sub = "../../6 Growth Models/Volume/umat-cauchy11.f"
# sub_start = 67

##### Orthotropic growth (only in 11 direction) with Bayly's evolution function #########
# sdvini = "sdvini.f"
# umat = "callumat_ortho.f"
# umat_cort = "../../6 Growth Models/Morphogenetic/Fiber/umat_length_linear.f"
# cor_start = 36
# umat_intro = "brain-ortho.f"
# umat_sub = "../../6 Growth Models/Orthotropic/umat_ortho_Bayly.f"
# sub_start = 77

##### Orthotropic growth (only in 11 direction) with linear evolution function #########
# sdvini = "sdvini.f"
# umat = "callumat_ortho.f"
# umat_cort = "../../6 Growth Models/Morphogenetic/Fiber/umat_length_linear.f"
# cor_start = 36
# umat_intro = "brain-ortho.f"
# umat_sub = "../../6 Growth Models/Orthotropic/umat_ortho_linear.f"
# sub_start = 77

##### Fiber growth with Bayly's evolution function ########################
# sdvini = "sdvini.f"
# umat = "callumat_fiber.f"
# umat_cort = "../../6 Growth Models/Morphogenetic/Fiber/umat_length_linear.f"
# cor_start = 36
# umat_intro = "brain-fiber.f"
# umat_sub = "../../6 Growth Models/Fiber/umat_fiber_Bayly.f"
# sub_start = 71

##### Fiber growth with Bayly's evolution function (pos or neg) ############
# sdvini = "sdvini.f"
# umat = "callumat_fiber.f"
# umat_cort = "../../6 Growth Models/Morphogenetic/Fiber/umat_length_linear.f"
# cor_start = 36
# umat_intro = "brain-fiber-neg.f"
# umat_sub = "../../6 Growth Models/Fiber/umat_fiber_Bayly_neg.f"
# sub_start = 72

##### Fiber growth with linear evolution function ########################
# sdvini = "sdvini.f"
# umat = "callumat_fiber.f"
# umat_cort = "../../6 Growth Models/Morphogenetic/Fiber/umat_length_linear.f"
# cor_start = 36
# umat_intro = "brain-fiber.f"
# umat_sub = "../../6 Growth Models/Fiber/umat_fiber_linear.f"
# sub_start = 70

##### Fiber growth with linear evolution function (pos or neg) ############
sdvini = "sdvini.f"
umat = "callumat_fiber.f"
umat_cort = "../../6 Growth Models/Morphogenetic/Fiber/umat_length_linear.f"
cor_start = 36
umat_intro = "brain-fiber-neg-linear.f"
umat_sub = "../../6 Growth Models/Fiber/umat_fiber_linear_neg.f"
sub_start = 70


##### Orientation #######################################################
orient = "../../6 Growth Models/Orient/orient-axons-r" + str(orientation) + ".f"

##### Perturbation #######################################################
# perturb[0] = '				if ((abs(coords(1)-center).lt.delta).and.(time(2).lt.tf)) then \n'
# perturb[1] = '				     R(i) = R(i) - perturb*time(2) \n'
# perturb[2] = '				endif \n\n'

##### UMAT ############################################################
f = open(umat_filename,"w")

# sdvini subroutine
a = open(sdvini,"r")
f.write(a.read())

# umat subroutine
b = open(umat,"r")
f.write(b.read())

# length morphogenetic growth for cortex
f.write("c...  ------------------------------------------------------------------\n")
f.write("      subroutine umat_cortex(stress,statev,ddsdde,sse,time,dtime,props,dfgrd1,ntens,ndi,nstatv,nprops,noel,npt,kstep,kinc) \n")
c = open(umat_cort,"r")
lines = c.readlines()
for line in lines[cor_start:]:
	f.write(line)

# umat for subcortex
f.write("\n\n\n")
r = open(umat_intro,"r")
f.write(r.read())

s = open(umat_sub,"r")
lines = s.readlines()
for line in lines[sub_start:]:
	if len(line) > 2 and line[-2] == "\r":
			line = line[0:-2] + "\n"
	# perturbation 
# 	if 'theg(i) = theg(i) - R(i) / dR(i)' in line:
# 		f.write(perturb[0])
# 		f.write(perturb[1])
# 		f.write(perturb[2])
	f.write(line)
	
# orient subroutine
o = open(orient,"r")
lines = o.readlines()
for line in lines:
	if len(line) > 2 and line[-2] == "\r":
			line = line[0:-2] + "\n"
	f.write(line)

f.close()