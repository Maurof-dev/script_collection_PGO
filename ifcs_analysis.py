import numpy as np
import pylab as plt
from math import sqrt

# The scope of the script is to plot the second order self-force constants (∂F_i/∂x_i, ∂F_i/∂y_i, ∂F_i/∂z_i) as a function of the atom --> to detect instabilities, assuming the self constribution is the most important

# The input file is the ifc.out output from ANADDB. It is assumed that the FULL force constant matrix is printed to file



def diagonal(namefile):
    file = open(namefile)
    data = file.readlines()
    
    atom_label = np.array([data[5*i*58].split() for i in range(0,57)])
    atom_label = np.array([float(data[5*i*58].split()[0]) for i in range(0,57)])

    IFCS_self = [data[5*j*58:5*j*58+5] for j in range(0,57)]

    
    IFCS_x = np.array([float(IFCS_self[i][2].split()[0]) for i in range(0,57)])
    IFCS_y = np.array([float(IFCS_self[i][3].split()[1]) for i in range(0,57)])
    IFCS_z = np.array([float(IFCS_self[i][4].split()[2]) for i in range(0,57)])

    IFCS_trace = np.array([float(IFCS_self[i][2].split()[0])+float(IFCS_self[i][3].split()[1])+float(IFCS_self[i][4].split()[2]) for i in range(0,57)])
    return [atom_label,IFCS_x,IFCS_y,IFCS_z,IFCS_trace]




    
def first_ifcs_line(namefile,centre):
    file = open(namefile)
    data = file.readlines()
    
    atom_label = np.array([data[5*i*58-5*i].split() for i in range(0,57)])
    atom_label = np.array([float(data[5*i*58 -5*i].split()[0]) for i in range(0,57)])

    centre = 57*centre 
    IFCS_self = [data[5*j:5*j+5] for j in range(0+centre,57+centre)]
    new_label = np.array([[float(IFCS_self[j][0].split()[0]),float(IFCS_self[j][0].split()[1])] for j in range(0,57)])
    pos_x = np.array([float(IFCS_self[i][1].split()[0]) for i in range(0,57)])
    pos_y = np.array([float(IFCS_self[i][1].split()[1]) for i in range(0,57)])
    pos_z = np.array([float(IFCS_self[i][1].split()[2]) for i in range(0,57)])
    pos_x0 = pos_x[0]
    pos_y0 = pos_y[0]
    pos_z0 = pos_z[0]

    pos_x = pos_x - pos_x0
    pos_y = pos_y - pos_y0
    pos_z = pos_z - pos_z0

    dist_list = np.array([sqrt(pos_x[i]**2 + pos_y[i]**2 + pos_z[i]**2) for i in range(0,57)]) # distance of all the 57 ions from the (1,1) reference atom
    
    IFCS_x = np.array([float(IFCS_self[i][2].split()[0]) for i in range(0,57)])
    IFCS_y = np.array([float(IFCS_self[i][3].split()[1]) for i in range(0,57)])
    IFCS_z = np.array([float(IFCS_self[i][4].split()[2]) for i in range(0,57)])

    IFCS_trace = np.array([float(IFCS_self[i][2].split()[0])+float(IFCS_self[i][3].split()[1])+float(IFCS_self[i][4].split()[2]) for i in range(0,57)])
    return [atom_label,IFCS_x,IFCS_y,IFCS_z,dist_list,IFCS_trace]#,IFCS_x,IFCS_y,IFCS_z,IFCS_trace]


    
def force_constants(namefile,flag,centre):
    if flag == 0:
        y = diagonal(namefile)
    if flag == 1:
        y = first_ifcs_line(namefile,centre)
    return y


PE_file = 'PGO_PE_ifcs.out'
FE_file = 'PGO_FErel_ifcs.out'
datafile = PE_file

check0 = force_constants(datafile,1,0)
check1 = force_constants(datafile,1,6)
check2 = force_constants(datafile,1,9)
check3 = force_constants(datafile,1,12)

check4 = force_constants(datafile,1,18)
check5 = force_constants(datafile,1,19)
check6 = force_constants(datafile,1,21)
check7 = force_constants(datafile,1,22)



check8 = force_constants(datafile,1,24)
check9 = force_constants(datafile,1,30)
check10 = force_constants(datafile,1,33)
check11 = force_constants(datafile,1,36)

check12 = force_constants(datafile,1,42)
check13 = force_constants(datafile,1,48)
check14 = force_constants(datafile,1,51)


'''
plt.plot(check[4],check[1],'-r.',label='xx')
plt.plot(check[4],check[2],'-b.', label='yy')
plt.plot(check[4],check[3],'-g.',label='zz')


plt.plot(check1[0],check1[1],'--r*')#,label='along x, FE')
plt.plot(check1[0],check1[2],'--b*')#, label='along y, FE')
plt.plot(check1[0],check1[3],'--g*')#,label='along z, FE')

plt.plot(check[0],check[4],'*r',label='trace')
plt.plot(check1[0],check1[4],'.c',label='trace, FE')


plt.plot(check[0],zeroline,'-m')

plt.legend()
plt.xlabel('atom index')
plt.ylabel('IFCs')
plt.show()
'''

plt.plot(check0[4],check0[5],'--g*',label='Ge(6l)')
plt.plot(check1[4],check1[5],'--b*',label='Ge(3k)')
'''
plt.plot(check2[4],check2[5],'--y*',label='Pb(3k)')
plt.plot(check3[4],check3[5],'--r*',label='Pb(6l)')
plt.plot(check4[4],check4[5],'--k*',label='Pb(1e)')

plt.plot(check5[4],check5[5],'-g*',label='Pb(2i)')
plt.plot(check6[4],check6[5],'-m*',label='Pb(1c)')
plt.plot(check7[4],check7[5],'-c*',label='Pb(2h)')
'''

plt.plot(check8[4],check8[5],'-m*',label='O(6l)')
plt.plot(check9[4],check9[5],'-c*',label='O(3k)')
plt.plot(check10[4],check10[5],'-y*',label='O(3k)')
plt.plot(check11[4],check11[5],'-r*',label='O(6l)')
plt.plot(check12[4],check12[5],'-k*',label='O(6l)')


plt.plot(check13[4],check13[5],'--m*',label='O(3j)')
plt.plot(check14[4],check14[5],'--c*',label='O(6l)')




plt.legend()
plt.ylabel('Tr[IFCs]')
plt.xlabel('distance from centre (Bohr)')
plt.show()





data = {'6l':check0[5][0],'3k':check1[5][0],'3k':check2[5][0],'6l':check3[5][0],'1e':check4[5][0],'2i':check5[5][0],'1c':check6[5][0],'2h':check7[5][0],'6l':check8[5][0],'3k':check9[5][0],'3k':check10[5][0],'6l':check11[5][0],'6l':check12[5][0],'3j':check13[5][0],'6l':check14[5][0]}


WPs = ['Ge(6l)','Ge(3k)', 'Pb(3k)','Pb(6l)','Pb(1e)','Pb(2i)','Pb(1c)','Pb(2h)','O1(6l)','O1(3k)','O2(3k)','O2(6l)','O3(6l)','O(3j)','O4(6l)']
IFCs = [check0[5][0],check1[5][0],check2[5][0],check3[5][0],check4[5][0],check5[5][0],check6[5][0],check7[5][0],check8[5][0],check9[5][0],check10[5][0],check11[5][0],check12[5][0],check13[5][0],check14[5][0]]

WPs_Ge = ['6l','3k']
WPs_Pb = ['3k','6l','1e','2i','1c','2h']
WPs_O = ['6l-1','3k-1','3k-2','6l-2','6l-3','3j','6l-4']


IFCs_Ge = [check0[5][0],check1[5][0]]
IFCs_Pb = [check2[5][0],check3[5][0],check4[5][0],check5[5][0],check6[5][0],check7[5][0]]
IFCs_O = [check8[5][0],check9[5][0],check10[5][0],check11[5][0],check12[5][0],check13[5][0],check14[5][0]]

bar_labels = ['Ge','Pb','O']
             
bar_colours = ['tab:red','tab:red','tab:blue','tab:blue','tab:blue','tab:blue','tab:blue','tab:blue','tab:green','tab:green','tab:green','tab:green','tab:green','tab:green','tab:green']
    
courses = list(data.keys())
values = list(data.values())

'''
fig = plt.figure(figsize = (9, 5))
 
# creating the bar plot
plt.bar(courses, values, color ='orange',
        width = 0.1)


 
#plt.xlabel("Courses offered")
plt.ylabel('on-site Tr[IFCs]')
plt.title("on-site interatomic force constants")
plt.show()
'''

fig, ax = plt.subplots()
ax.bar(WPs_Ge, IFCs_Ge, label='Ge', color='red')
ax.bar(WPs_Pb, IFCs_Pb, label='Pb', color='orange')
ax.bar(WPs_O, IFCs_O, label='O', color='blue')

ax.set_ylabel('on-site Tr[IFCs]')
ax.set_title('on-site interatomic force constants')
ax.legend()

plt.show()
