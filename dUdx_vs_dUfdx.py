# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# File name: dUdx_vs_dUfdx.py
# Created by: gemusia
# Creation date: 13-09-2017
# Last modified: 14-09-2017 14:24:26
# Purpose: computes and creates plot
#          comparing terms
#          (V-U)dU/dx with (V-U)dUf/dx
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




import apriori_tests_timestat as at


def f_dUdx_vs_dUfdx(V1,V2,V3,U1,U2,U3,dUdx,dUdy,dUdz,dUfdx,dUfdy,dUfdz):
    DNS_derivative =  (V1-U1)*dUdx + (V2-U2)*dUdy + (V3-U3)*dUdz
    LES_derivative =  (V1-U1)*dUfdx + (V2-U2)*dUfdy + (V3-U3)*dUfdz
    return (DNS_derivative - LES_derivative)/DNS_derivative


dUdx_vs_dUfdx = [['Vx','Vy','Vz','Ux','Uy','Uz','dUxdx','dUxdy','dUxdz','dUfxdx','dUfxdy','dUfxdz'],
                 ['Vx','Vy','Vz','Ux','Uy','Uz','dUydx','dUydy','dUydz','dUfydx','dUfydy','dUfydz'],
                 ['Vx','Vy','Vz','Ux','Uy','Uz','dUzdx','dUzdy','dUzdz','dUfzdx','dUfzdy','dUfzdz']]

if __name__=='__main__':
    for StNo in at.ptype:

        print "ptype = ",at.ptype
        for stattype in ("pmean","pstd"):

            # STATISTICS
            pstat_test3 = at.pfields.equationP(StNo,f_dUdx_vs_dUfdx,stattype,"symm",*dUdx_vs_dUfdx)  


            # velocity statistics


            for pKey_test3 in at.keys_no_yplus(pstat_test3.keys()):

                ptermfullvel = at.hfig.Homfig(title="pterm ", ylabel="$((V-U)_j*dU/dx_j)^{2}$")
                plotFileNamePterm = at.pict_path + "dUdx-dUfdx_"+stattype+ "_"+StNo+".eps"

                for direction in range(3):
                    ptermfullvel.add_plot(pstat_test3["yplus"],pstat_test3[pKey_test3]/at.termplus,linestyle='dotted',label='% err($(V-U)_j*d\overline{U}/dx_j$')
                    
                ptermfullvel.hdraw()
                ptermfullvel.save(plotFileNamePterm)
                print "plot created: " + plotFileNamePterm
                at.plt.close(ptermfullvel.fig)
