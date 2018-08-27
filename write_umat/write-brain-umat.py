import sys

##### PRELIM ##########################################################

name = sys.argv[1]
orientation = sys.argv[2]
umat_filename = "../umat_" + name + ".f"

##### Fiber growth with linear evolution function (pos or neg) ############
sdvini = "sdvini.f"
umat = "callumat_fiber.f"
umat_cort = "../../6 Growth Models/Morphogenetic/Fiber/umat_length_linearmorph.f"
cor_start = 36
umat_intro = "brain-fiber-neg-linear.f"
orient = "../../6 Growth Models/Orient/axon" + str(orientation) + ".f"
umat_sub = "../../6 Growth Models/Fiber/umat_subcortex.f"
sub_start = 76


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

# umat heading for subcortex
f.write("\n\n\n")
r = open(umat_intro,"r")
f.write(r.read())

# orient subroutine
o = open(orient,"r")
f.write(o.read())
f.write('\n')

# umat body for subcortex
s = open(umat_sub,"r")
lines = s.readlines()
for line in lines[sub_start:]:
   if len(line) > 2 and line[-2] == "\r":
         line = line[0:-2] + "\n"
   f.write(line)
   


f.close()