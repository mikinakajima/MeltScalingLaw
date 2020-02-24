import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys
import matplotlib.colors as colors
import matplotlib.cm as cmx

#I get soome errors but this still works...
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def test_figure3(xlabel,ylabel,title,x,y,ynum,xlim_in,ylim_in,xint,yint,ydata_label,font_label,font_legend,font_xylabels,border_width,alpha_num,line_width, fig, ax, form, line_type,vesc, *args):


    #Cvalues, scalarMap = color_map(3, 'hot') #tab20c

    #print('a')
    print(len(args))

    #plt.style.use('seaborn-notebook')

    text=[]
    mantle_high_v=0
    
    if(len(args)>0):
        xfit=args[0] #fitting
        yfit=args[1] #fitting

    if(len(args)>2):
        text=args[2]
        text_loc=args[3]
        text_font=args[4]
    if(len(args)>5):
        mantle_high_v=args[5]
        print(mantle_high_v,'hi')


    
    ax.set_color_cycle(sns.color_palette("coolwarm",len(x)))    
    ax.set(xlim=xlim_in,ylim=ylim_in)
    ax.set_xlabel(xlabel,fontsize=font_label)
    ax.set_ylabel(ylabel,fontsize=font_label)
    ax.set_title(title,fontsize=font_label)
    
    
    ax.xaxis.set_ticks(np.arange(xlim_in[0],xlim_in[1]+0.01, xint))
    ax.yaxis.set_ticks(np.arange(ylim_in[0],ylim_in[1]+0.01, yint))

    ax.xaxis.set_tick_params(labelsize=font_xylabels)
    ax.yaxis.set_tick_params(labelsize=font_xylabels)

    print(x.shape,'b')
    
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(border_width)

    #ax.set_color_cycle(sns.color_palette("coolwarm", len(x)))

    if form=="loglog":
        for i in range(0,len(x)):
            ax.loglog(x[i][:],y[i][:],line_type[i],label=ydata_label[i], linewidth=line_width[i],alpha=alpha_num)
            ax.legend(fontsize=font_legend,frameon=False)
    elif form=="semilogx":
        #print("hello")
        ax.semilogx(x,y,label=ydata_label, linewidth=line_width,alpha=alpha_num)
    else:
        if len(ydata_label)==0:

            if vesc<0.5:
            
                gamma_01=[0, 1, 2, 3, 6, 7, 10, 11, 12, 13, 14]
                gamma_05=[4]
                gamma_003=[5]
                gamma_02=[8]
            
                ax.set_color_cycle(sns.color_palette("summer",len(gamma_01)))
                for j in range(0, len(gamma_01)):
                    i=gamma_01[j]
                    print(i,'a')
                    ax.plot(x[i][:],y[i][:],line_type[i])

                ax.set_color_cycle(sns.color_palette("pink",len(gamma_05)))
                for j in range(0, len(gamma_05)):
                    i=gamma_05[j]
                    ax.plot(x[i][:],y[i][:],line_type[i])

                ax.set_color_cycle(sns.color_palette("ocean",len(gamma_003)))
                for j in range(0, len(gamma_003)):
                    i=gamma_003[j]
                    ax.plot(x[i][:],y[i][:],line_type[i])

                ax.set_color_cycle(sns.color_palette("PuRd",len(gamma_02)))
                for j in range(0, len(gamma_02)):
                    i=gamma_02[j]
                    ax.plot(x[i][:],y[i][:],line_type[i])
            else:
                gamma_01=[1, 2, 3, 4, 5, 6]
                gamma_05=[0]
                gamma_03=[7, 8, 9, 10]

                ax.set_color_cycle(sns.color_palette("summer",len(gamma_01)))
                for j in range(0, len(gamma_01)):
                    i=gamma_01[j]
                    print(i,'c')
                    ax.plot(x[i][:],y[i][:],line_type[i])

                ax.set_color_cycle(sns.color_palette("pink",len(gamma_05)))
                for j in range(0, len(gamma_05)):
                    i=gamma_05[j]
                    ax.plot(x[i][:],y[i][:],line_type[i])

                ax.set_color_cycle(sns.color_palette("RdPu",len(gamma_03)))
                for j in range(0, len(gamma_03)):
                    i=gamma_03[j]
                    ax.plot(x[i][:],y[i][:],line_type[i])
                                
        else:

            if vesc<0.5:
            
                gamma_01=[0, 1, 2, 3, 6, 7, 10, 11, 12, 13, 14]
                gamma_05=[4]
                gamma_003=[5]
                gamma_02=[8]
            
                ax.set_color_cycle(sns.color_palette("summer",len(gamma_01)))
                for j in range(0, len(gamma_01)):
                    i=gamma_01[j]
                    print(i,'a')
                    ax.plot(x[i][:],y[i][:],line_type[i], label=ydata_label[i])

                ax.set_color_cycle(sns.color_palette("pink",len(gamma_05)))
                for j in range(0, len(gamma_05)):
                    i=gamma_05[j]
                    ax.plot(x[i][:],y[i][:],line_type[i], label=ydata_label[i])

                ax.set_color_cycle(sns.color_palette("ocean",len(gamma_003)))
                for j in range(0, len(gamma_003)):
                    i=gamma_003[j]
                    ax.plot(x[i][:],y[i][:],line_type[i],label=ydata_label[i])

                ax.set_color_cycle(sns.color_palette("PuRd",len(gamma_02)))
                for j in range(0, len(gamma_02)):
                    i=gamma_02[j]
                    ax.plot(x[i][:],y[i][:],line_type[i],label=ydata_label[i])
            else:
                gamma_01=[1, 2, 3, 4, 5, 6]
                gamma_05=[0]
                gamma_03=[7, 8, 9, 10]

                ax.set_color_cycle(sns.color_palette("summer",len(gamma_01)))
                for j in range(0, len(gamma_01)):
                    i=gamma_01[j]
                    print(i,'c')
                    ax.plot(x[i][:],y[i][:],line_type[i], label=ydata_label[i])

                ax.set_color_cycle(sns.color_palette("pink",len(gamma_05)))
                for j in range(0, len(gamma_05)):
                    i=gamma_05[j]
                    ax.plot(x[i][:],y[i][:],line_type[i], label=ydata_label[i])

                ax.set_color_cycle(sns.color_palette("RdPu",len(gamma_03)))
                for j in range(0, len(gamma_03)):
                    i=gamma_03[j]
                    ax.plot(x[i][:],y[i][:],line_type[i], label=ydata_label[i])

            

            
            #for i in range(0,len(x)):
            #    ax.plot(x[i][:],y[i][:],line_type[i],label=ydata_label[i])

        if len(args)==0:
            pass
            #print("hi")
        elif mantle_high_v==0:
            ax.plot(xfit,yfit,line_type[i],linewidth=3,color='black')
        else:
            ax.set_color_cycle(sns.color_palette("coolwarm",len(x)))
            #sns.palplot(sns.color_palette("Blues"))

            
            for i in range(0,len(x)): #len(x)
                ax.plot(xfit[i][:],yfit[i][:],'--', linewidth=line_width,alpha=alpha_num)

                

    ax.legend(bbox_to_anchor=(0., -0.1, 1., -.2), #loc=3,
                  ncol=4,  mode="expand", borderaxespad=0., frameon=False, fontsize=font_legend)


        #ax.legend(fontsize=font_legend,frameon=False)
    
    if len(text)>0:
        ax.text(xlim_in[0]+text_loc[0]*(xlim_in[1]-xlim_in[0]),ylim_in[0]+text_loc[1]*(ylim_in[1]-ylim_in[0]),text,fontsize=text_font)



    



    #print(x)
    #print(y)



