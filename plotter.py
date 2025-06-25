from IPython.core.pylabtools import figsize
from event import Event
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import matplotlib.patches as mpatches
from numpy.ma.core import asarray


class Plotter:

    def __init__(self, event):
        print("plotter: initilization")
        self.event = event
        #self.Jump(0, 0.5)

    def Collect(self, dQthreshold):
        self.xx = np.array([])
        self.yy = np.array([])
        self.zz = np.array([])
        self.cc = np.array([])
        self.tt = np.array([])
        self.ee = np.array([])
        self.ll = np.array([])
        self.aa = np.array([])
        self.ii = np.array([])
        self.pp = np.array([])
        self.pp_2 = np.array([])
        self.USER_COLORS = ['black', 'red', 'blue', 'magenta']

        nColor = len(self.USER_COLORS)
        mm2m = 0.001
        mm2cm = 0.1
        for i, track in enumerate(self.event.tracks):
            depoList = track.association['depoList']
            raw_muon = self.event.read_muon_direction()
            if raw_muon is None or np.all(raw_muon == 0):
                raise RuntimeError("no muon direction!")
            muon_direction = raw_muon / np.linalg.norm(raw_muon)
            ancestor = track.association['ancestor']
            ancestor = self.event.tracks[ancestor].GetPDGCode()
            for di in depoList:
                depo = self.event.depos[di]
                x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                t = (depo.GetStart().T() + depo.GetStop().T()) / 2 # ns
                e = depo.GetEnergyDeposit() # MeV
                l = depo.GetTrackLength() *mm2cm # most are 0.5 cm
                edep_vec = np.array([x, y, z])
                cos = np.dot(edep_vec, muon_direction) / (np.linalg.norm(edep_vec) * np.linalg.norm(muon_direction))
                aa = abs(np.degrees(np.arccos(np.clip(cos, -1, 1))))
                if ancestor==2112:
                    ii=0
                elif ancestor in (13,-13):
                    ii=1
                else:
                    ii=2
                if ii==0 and aa > 10:
                    ##which means a neutron edep we are able to find
                    pp=0
                elif ii==0 and aa<=10:
                    pp=1
                    ##which means a neutron that is wrongly excluded
                elif ii==1 and aa > 10:
                    pp=2
                    ##which means a muon wrongly included
                elif ii==1 and aa<=10:
                    pp=3
                    ##which means a muon successfully excluded
                elif ii==2 and aa > 10:
                    pp=4
                    ##which means another particle wrongly included
                elif ii==2 and aa<=10:
                    pp=5
                else:
                    pp=6
                print(f"pp is {pp}")


                if ii==0 and aa > 20:
                    ##which means a neutron edep we are able to find
                    pp_2=0
                elif ii==0 and aa<=20:
                    pp_2=1
                    ##which means a neutron that is wrongly excluded
                elif ii==1 and aa > 20:
                    pp_2=2
                    ##which means a muon wrongly included
                elif ii==1 and aa<=20:
                    pp_2=3
                    ##which means a muon successfully excluded
                elif ii==2 and aa > 20:
                    pp_2=4
                    ##which means another particle wrongly included
                elif ii==2 and aa<=20:
                    pp_2=5
                else:
                    pp_2=6
                print(f"pp_2 is {pp}")
                if self.event.ChargeBirksLaw(e, l) > dQthreshold: # detector threshold
                    self.xx = np.append(self.xx, x)
                    self.yy = np.append(self.yy, y)
                    self.zz = np.append(self.zz, z)
                    self.tt = np.append(self.tt, t)
                    self.ee = np.append(self.ee, e)
                    self.ll = np.append(self.ll, l)
                    self.cc = np.append(self.cc, self.USER_COLORS[ancestor % nColor])
                    self.aa =np.append(self.aa, aa)
                    self.ii = np.append(self.ii, ii)
                    self.pp = np.append(self.pp, pp)
                    self.pp_2 = np.append(self.pp_2, pp_2)
    def Jump(self, entryNo, dQthreshold):
        self.event.Jump(entryNo)
        self.Collect(dQthreshold)

    def Next(self, dQthreshold):
        self.event.Next()
        self.Collect(dQthreshold)

    def Prev(self, dQthreshold):
        self.event.Prev()
        self.Collect(dQthreshold)

    def Draw(self, axis='yz', value='time', energy='GeV',
             markerSize=3, cmap='jet', vmax=2000):
        # load event & collect deposit info
        self.event.Jump(self.event.currentEntry)
        self.Collect(0)
        mapping = {'x': self.xx, 'y': self.yy, 'z': self.zz, 'e': self.ee, 't': self.tt}
        mask = self.ee >= 0.5

        # 1) create 3×2 grid, then delete the last two axes (positions 5 and 6)
        fig, axes = plt.subplots(3, 2, figsize=(10, 12), dpi=100)
        axes = axes.flatten()


        # 2) reassign only the first four axes
        ax1, ax2, ax3, ax4, ax5 ,ax6 = axes[:6]

        fig.suptitle(self.event.vertex.GetReaction())



        # — Panel 1: time scatter —
        plot_time = ax1.scatter(
            mapping[axis[1]], mapping[axis[0]],
            c=self.tt,
            cmap='viridis',
            norm=matplotlib.colors.LogNorm(),
            s=15,
            edgecolors='none',
            linewidths=1,
            zorder=1
        )

        ax1.set_xlabel('Time (ns)')
        ax1.set_title('Time [ns]')
        fig.colorbar(plot_time, orientation='vertical',label="time[ns]",ax=ax1,pad=0.02,shrink=1.0)

        # — Panel 2: energy scatter —
        plot_energy = ax2.scatter(
            mapping[axis[1]], mapping[axis[0]],
            c=self.ee * 2,
            cmap='viridis',
            vmax=2.5,
            s=15,
            edgecolors='none',
            linewidths=1,
            zorder=1
        )
        x_start, x_end = ax2.get_xlim()
        y_start, y_end = ax2.get_ylim()

        ax2.set_xlabel('Energy [GeV]')
        ax2.set_title('Energy [MeV/cm]')
        fig.colorbar(plot_energy, orientation='vertical',label="MeV/cm", ax=ax2,pad=0.02,shrink=1.0)


        # — Panel 3: neutron direction + filtered scatter —
        dir_vec = self.event.read_neutron_direction()
        if dir_vec is None:
            print(f"Event {self.event.currentEntry}: no neutron direction, skipping arrow.")
        else:
            comp = {'x': 0, 'y': 1, 'z': 2}
            dx = dir_vec[comp[axis[1]]]
            dy = dir_vec[comp[axis[0]]]
            vec = np.array([dx, dy])
            vec = vec / np.linalg.norm(vec)
            ax3.arrow(
                0, 0,
                vec[0], vec[1],
                color='blue',
                length_includes_head=True,
                head_width=0.1,
                head_length=0.1
            )
            ax4.arrow(
                0, 0,
                vec[0], vec[1],
                color='blue',
                length_includes_head=True,
                head_width=0.1,
                head_length=0.1
            )

        ax3.scatter(
            mapping[axis[1]][mask],
            mapping[axis[0]][mask],
            c=self.cc[mask],
            s=markerSize
        )
        ax3.set_xlabel('Ancestor ID(y-z plane, [m])')

        ax4.scatter(
            mapping[axis[1]][mask],
            mapping[axis[0]][mask],
            c=self.cc[mask],
            s=markerSize
        )
        ax4.set_xlabel('Ancestor ID')
        ax4.set_xlim(x_start,x_end)
        ax4.set_ylim(y_start,y_end)
        ax4.set_xlabel("zoom out version of selected edep(y-z plane, [m])")




        data = np.asarray(self.event.edep_based_information())
        if data is None or len(data) == 0:
            print("no data, skipping")
            return
        arr = np.asarray(data)
        if  any(len(row) < 3 for row in data):
            print(f"Event {self.event.currentEntry}:imcomplete information, skip and return")
            return

        all_length = data[:, 0]
        all_energy = data[:, 1]
        all_angle = data[:, 2]
        plot_energy_angle = ax5.scatter(
            all_energy, all_angle,
            c=all_length,
            cmap='viridis',
            s=15,
            edgecolors='none',
            linewidths=1,
            zorder=1
        )
        ax5.set_xlim(0,3)
        ax5.set_ylim(0,3.1416)
        fig.colorbar(plot_energy_angle, label="length[cm]", ax=ax5,pad=0.02,shrink=1.0)
        ax5.set_xlabel("energy[MeV]")
        ax5.set_ylabel("angle[rad]")


        plot_length_angle = ax6.scatter(
            all_length, all_angle,
            c=all_energy,
            cmap='viridis',
            s=15,
            edgecolors='none',
            linewidths=1,
            zorder=1
        )
        ax6.set_xlim(0,400)
        ax6.set_ylim(0,3.1416)
        fig.colorbar(plot_length_angle, label="energy[MeV]", ax=ax6,pad=0.02,shrink=1.0)
        ax6.set_xlabel("Length[cm]")
        ax6.set_ylabel("angle[rad]")



        # Rest of your styling, text labels, and save logic unchanged:
        for ax in (ax1, ax2):
            if self.event.evgen == 'Genie':
                ax.set_ylim(-4, 4)
                ax.set_xlim(-2, 8)
                xpos, ypos, interval = -1.8, 3.6, 0.35
            elif self.event.evgen == 'Marley':
                ax.set_ylim(-1, 1)
                ax.set_xlim(-1, 0.5)
                xpos, ypos, interval = -0.95, 0.9, 0.1
            else:
                print("Unknown event generator!")
                sys.exit()
            ax.tick_params(axis='both', direction='in', length=2)
            ax.set_xlabel(f'{axis[1]} [m]')
            if ax is ax1:
                ax.set_ylabel(f'{axis[0]} [m]')

        nColor, countnegId = len(self.USER_COLORS), 0
        entries = []
        for i, particle in enumerate(self.event.vertex.Particles):
            trkId = particle.GetTrackId()
            if trkId < 0:
                countnegId += 1
                continue
            mom = particle.GetMomentum()
            KE = mom.E() - mom.M()
            name = f'{particle.GetName()}: {KE :.1f} MeV'
            color = self.USER_COLORS[(i - countnegId) % nColor]
            # now annotate on ax3 in axes‐fraction coords
            x_frac = 0.02  # 2% in from the left
            y_start = 0.98  # start 98% up from bottom
            y_gap_frac = 0.05  # 5% vertical spacing
            entries.append((name, color))

            for idx,(text,color) in enumerate(entries):
                ax3.text(
                    x_frac,
                    y_start - idx * y_gap_frac,
                    text,
                    transform=ax3.transAxes,
                    va='top',
                    color=color,
                    fontsize=8
                )
        fig.savefig(
            f"{self.event.plotpath}/particle_{value}_{axis}_{energy}_evt_{self.event.currentEntry}.pdf"
        )
        plt.clf()
        plt.close()

    # For GeV events
    def hist_dEdx(self):
        plt.hist(self.ee*2, range=(0,6), bins=100)
        plt.xlabel('dE/dx [MeV/cm]')
        plt.draw()
        plt.savefig(self.event.plotpath + '/dE_dx_evt_%d.pdf' % self.event.currentEntry)
        plt.clf() # important to clear figure
        plt.close()
        #plt.show()

    # For low E MeV events, edep length is smaller than 0.5cm (5mm step limit in gdml)
    def hist_LowE_dEdx(self):
        plt.hist(np.divide(self.ee, self.ll), range=(0,30), bins=500)
        plt.xlabel('dE/dx [MeV/cm]')
        plt.draw()
        plt.savefig(self.event.plotpath + '/LowE_dE_dx_evt_%d.pdf' % self.event.currentEntry)
        plt.clf() # important to clear figure
        plt.close()



    def muon_edep_histogram(self):
        raw_angle=self.aa
        angle= np.array(raw_angle)
        muon_angle = []
        for i in range (len(angle)):
            if np.array(self.ii)[i]==1:
                muon_angle.append(angle[i])


        # fig , ax = plt.subplots(8,5,dpi=300)
        # counts =
        # mid_bar =
        data_min = np.min(muon_angle)
        data_max = np.max(muon_angle)
        x_min_plot = data_min - (data_max - data_min) * 0.05
        x_max_plot = data_max + (data_max - data_min) * 0.05
        bins_count=50
        bins = np.linspace(x_min_plot, x_max_plot, bins_count)
        plt.hist(muon_angle,bins=bins, density=False, color='skyblue', edgecolor='black',alpha=0.5)

        plt.xlim(x_min_plot, x_max_plot)
        plt.xlabel('edep relative angle to muon direction[degree]')
        plt.ylabel('count')
        plt.title("muon edep relative direction to muon's initial direction(0.2GeV antinumu) ")
        plt.savefig(self.event.plotpath + '/muon_edep_angle_histogram0.2GeV.pdf')

    def Draw_angle(self, axis='yz', value='charge', energy='GeV', markerSize=3, cmap='jet', vmax=2000):
            # particle, timing, dE/dx
        mapping = {'x': self.xx, 'y': self.yy, 'z': self.zz}

        fig, axes = plt.subplots(2, 3, figsize=(15, 8), dpi=600)
        ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()
        fig.suptitle(self.event.vertex.GetReaction())
        cb_ax = fig.add_axes([.94, .124, .02, .754])

            # particle plot
        ax1.scatter(mapping[axis[1]], mapping[axis[0]], c=self.cc, s=markerSize)
        if value == 'time':
            # timing plot
            norm = matplotlib.colors.LogNorm(vmin=1, vmax=vmax)
            plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.tt, cmap=cmap, vmin=1, vmax=vmax,norm=matplotlib.colors.LogNorm(), s=markerSize)
            cb_ax.set_xlabel('ns')

        elif value == 'charge':
            # charge plot
            norm = matplotlib.colors.LogNorm(vmin=1, vmax=vmax)
            if energy == 'GeV':
                plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.ee * 2, cmap=cmap, vmax=4,
                                          s=markerSize)
            elif energy == 'MeV':
                # plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=np.divide(self.ee, self.ll), cmap=cmap, vmax=4, s=markerSize)
                plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.ee, cmap=cmap, vmax=2.5,norm=norm,
                                          s=markerSize)
            cb_ax.set_xlabel('MeV')

        fig.colorbar(plot_12, orientation='vertical', cax=cb_ax)
            # fig.tight_layout()
        for ax in (ax1, ax2, ax3, ax4, ax5, ax6):
            if self.event.evgen == 'Genie':
                ax.set_ylim(-4, 4)
                ax.set_xlim(-2, 8)
                xpos = -1.8
                ypos = 3.6
                interval = 0.35
            else:
                print("Unknown event generator!")
                sys.exit()
            ax.tick_params(axis='y', direction='in', length=2)
            ax.tick_params(axis='x', direction='in', length=2)
            ax.set_xlabel(f'{axis[1]} [m]')
            if ax == ax1:
                ax.set_ylabel(f'{axis[0]} [m]')

        nColor = len(self.USER_COLORS)
        countnegId = 0
        for i, particle in enumerate(self.event.vertex.Particles):
                # Skip negative trk id: in the case of Marley events,
                # this usually is the final nucleus before deexcitation that G4 doesn't track
                # the kinematics are not correct either
            trkId = particle.GetTrackId()
            if trkId < 0:
                    # Because we skipped the trk id, need to subtract the index for the proper coloring of deposits
                countnegId += 1
                continue
            name = particle.GetName()
            color = self.USER_COLORS[(i - countnegId) % nColor]
                # pdg = particle.GetPDGCode()
                # name = particle.GetName()
                # trkId = particle.GetTrackId()
            mom = particle.GetMomentum()
            KE = mom.E() - mom.M()
            name = '%s: %.1f MeV' % (name, KE)

            ax1.text(xpos, ypos, name, color=color)
            ypos -= interval
        #####
        #draw the arrow
        comp = {'x': 0, 'y': 1, 'z': 2}

        dn = self.event.read_neutron_direction()
        dm = self.event.read_muon_direction()
        if dn is None or dm is None:
            return
        for vec, color, label in ((dn, 'blue', 'neutron'), (dm, 'red', 'muon')):
                v = np.array([vec[comp['z']], vec[comp['y']]])
                v = v / np.linalg.norm(v)
                ax2.arrow(0, 0, v[0], v[1],
                          color=color, length_includes_head=True,
                          head_width=0.1, head_length=0.1)
                # legend via patches
        handles = [
                mpatches.Patch(color='blue', label='neutron(arrow)'),
                mpatches.Patch(color='red', label='muon(arrow)'),
            ]
        ax2.legend(handles=handles, loc='upper right', fontsize=8)
        # ax2.plot()

        dathreshold=10
        mask=self.aa>dathreshold
        all_colors = np.where(self.ii == 0, 'red',np.where(self.ii == 1, 'black', 'grey'))
        print("self.ii is ",self.ii)
        ax3.scatter(
                mapping[axis[1]][mask],
                mapping[axis[0]][mask],
                s=markerSize,
                c=all_colors[mask],  # pick only the masked ones
        )
        handles_2 = [
                mpatches.Patch(color='red', label='neutron'),
                mpatches.Patch(color='black', label='misidentify muon'),
                mpatches.Patch(color='grey', label='other primary particles'),
        ]
        ax3.legend(handles=handles_2, loc='upper right', fontsize=8)
        ax3.set_title("particle survived under the filter(10 degrees)")

        conditions = [
                (self.pp == 0),
                (self.pp == 1),
                (self.pp == 2),

        ]
        choices = [
                'red',  # for 0
                'black',  # for 1
                'green',  # for 3

        ]
        conditions_6 = [
            (self.pp_2 == 0),
            (self.pp_2 == 1),
            (self.pp_2 == 2),
        ]

            # any element not matching one of the above gets 'grey'
        all_colors_3 = np.select(conditions, choices, default='grey')
        all_colors_6 = np.select(conditions_6,choices, default='grey')

        dathreshold_2=20
        mask2=self.aa>dathreshold_2
        ax4.scatter(mapping[axis[1]][mask2],
                    mapping[axis[0]][mask2],
                    s=markerSize,
                    c=all_colors[mask2])
        ax4.legend(handles=handles_2, loc='upper right', fontsize=8)
        ax4.set_title("particle survived by the filter(20 degrees)")
        mask3 = ~np.isin(self.pp, [3, 4, 5])
        mask6 = ~np.isin(self.pp_2, [3, 4, 5])

        ax5.scatter(mapping[axis[1]][mask3],
                    mapping[axis[0]][mask3],
                    s=markerSize,
                    c=all_colors_3[mask3]
                    )
        handles_3 = [
                mpatches.Patch(color='red', label='neutron edep able to find'),
                mpatches.Patch(color='black', label='neutron edep wrongly exclude'),
                mpatches.Patch(color='green', label='muon wrongly included'),
                mpatches.Patch(color='grey',label='other primary particle edep')
            ]
        ax5.legend(handles=handles_3, loc='upper right', fontsize=8)
        ax5.set_title("particle purity and efficiency(10 degrees)")

        ax6.scatter(mapping[axis[1]][mask6],
                        mapping[axis[0]][mask6],
                        s=markerSize,
                        c=all_colors_6[mask6]
                        )

        ax6.legend(handles=handles_3, loc='upper right', fontsize=8)
        ax6.set_title("particle purity and efficiency(20 degrees)")

        fig.savefig(self.event.plotpath + '/particle_%s_%s_%s_evt_%d.pdf' % (value, axis, energy, self.event.currentEntry))
        plt.clf()  # important to clear figure
        plt.close()
            # plt.show()

    def hist_dx(self):
        plt.hist(self.ll, range=(0,4), bins=100)
        plt.xlabel('dx [cm]')
        plt.draw()
        plt.savefig(self.event.plotpath + '/single_edep_dx_evt_%d.pdf' % self.event.currentEntry)
        plt.clf() # important to clear figure
        plt.close()
        return self.ll #cm

    def hist_trklength(self):
        all_tracks_length = np.array([])
        for i, track in enumerate(self.event.tracks):
            all_tracks_length = np.append(all_tracks_length, track.length['selfDepo'])

        plt.hist(all_tracks_length, range=(0,6), bins=100)
        plt.xlabel('track length [cm]')
        plt.draw()
        plt.savefig(self.event.plotpath + '/primary_track_length_evt_%d.pdf' % self.event.currentEntry)
        plt.clf() # important to clear figure
        plt.close()
        return all_tracks_length #cm

    def evt_maxdtdr(self): # per evt has one max dt
        nedeps = len(self.tt)
        max_edep_dr, max_edep_dt = 0, 0
        for iedep in range(nedeps):
            for jedep in range(iedep+1, nedeps):
                temp_dt = abs(np.float128(self.tt[iedep] - self.tt[jedep])) # ns
                # unit m
                temp_dr = math.sqrt( pow(np.float128(self.xx[iedep] - self.xx[jedep]), 2) + pow(np.float128(self.yy[iedep] - self.yy[jedep]), 2) + pow(np.float128(self.zz[iedep] - self.zz[jedep]), 2) )
                if temp_dt > max_edep_dt:
                    max_edep_dt = temp_dt
                if temp_dr > max_edep_dr:
                    max_edep_dr = temp_dr

        return [max_edep_dr*100, max_edep_dt] #cm, ns

    def plot_angle_distribution(self):
        # Original code preserved
        muon_distribution = []
        neutron_distribution = []
        relative_distribution = []
        n_events = min(10, self.event.nEntry)
        for idx in range(n_events):
            self.event.Jump(idx)
            neutron_direction = self.event.read_neutron_direction()
            muon_direction = self.event.read_muon_direction()
            if neutron_direction is None or muon_direction is None or np.all(neutron_direction)==0 or np.all(muon_direction)==0:
                print(f"Event {self.event.currentEntry}: no neutron direction, skipping arrow!!!.")
                continue
            neutron_theta = np.arctan((neutron_direction[0] ** 2 + neutron_direction[1] ** 2) ** 0.5 / neutron_direction[2])
            if neutron_theta<=0:
                neutron_theta=3.1415+neutron_theta
            if neutron_direction[0]>0:
                if neutron_direction[1]>0:
                    neutron_phi = np.arctan(neutron_direction[1] / neutron_direction[0])
                else:
                    neutron_phi = 6.2832+np.arctan(neutron_direction[1] / neutron_direction[0])
            elif neutron_direction[0]<0:
                if neutron_direction[1]>0:
                    neutron_phi = 3.14159+np.arctan(neutron_direction[1] / neutron_direction[0])
                else:
                    neutron_phi = 3.14159+np.arctan(neutron_direction[1] / neutron_direction[0])
            neutron_distribution.append([neutron_theta, neutron_phi])

            muon_theta = np.arctan((muon_direction[0] ** 2 + muon_direction[1] ** 2) ** 0.5 / muon_direction[2])
            muon_phi = 0
            if muon_theta <= 0:
                muon_theta = 3.1415 + muon_theta
            if muon_direction[0] > 0:
                if muon_direction[1] > 0:
                    muon_phi = np.arctan(muon_direction[1] / muon_direction[0])
                else:
                    muon_phi = 6.2832 + np.arctan(muon_direction[1] / muon_direction[0])
            elif muon_direction[0] < 0:
                if muon_direction[1] > 0:
                    muon_phi = 3.14159 + np.arctan(muon_direction[1] / muon_direction[0])
                else:
                    muon_phi = 3.14159 + np.arctan(muon_direction[1] / muon_direction[0])

            muon_distribution.append([muon_theta, muon_phi])
            relative_theta=muon_theta-neutron_theta
            relative_phi=muon_phi-neutron_phi
            if relative_theta<0:
                relative_theta=relative_theta+2*np.pi
            if relative_phi<0:
                relative_phi=relative_phi+2*np.pi
            colors = 'red' if self.event.selectmuonneutronevent(idx) == 1 else 'blue'
            relative_distribution.append((relative_theta, relative_phi, colors))

        # Convert to arrays and degrees
        neutron_thetas, neutron_phis = np.array(neutron_distribution).T
        muon_thetas, muon_phis = np.array(muon_distribution).T

        rel_thetas = [row[0] for row in relative_distribution]
        rel_phis = [row[1] for row in relative_distribution]
        colors_rel = [row[2] for row in relative_distribution]


        neutron_thetas = np.degrees(neutron_thetas)
        neutron_phis = np.degrees(neutron_phis)
        muon_thetas = np.degrees(muon_thetas)
        muon_phis = np.degrees(muon_phis)
        rel_thetas = np.degrees(rel_thetas)
        rel_phis = np.degrees(rel_phis)



        # Compute true separation angles
        separations = []
        for idx in range(n_events):
            self.event.Jump(idx)
            n_dir = self.event.read_neutron_direction()
            m_dir = self.event.read_muon_direction()
            if n_dir is None or m_dir is None or np.all(n_dir)==0 or np.all(m_dir)==0 :
                print(f"Event {self.event.currentEntry}: no neutron direction, skipping arrow.")
                continue
            cosang = np.dot(n_dir, m_dir) / (np.linalg.norm(n_dir) * np.linalg.norm(m_dir))
            cosang = np.clip(cosang, -1, 1)
            separations.append(np.degrees(np.arccos(cosang)))

        # Plotting
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))
        axs[0, 0].scatter(muon_thetas, muon_phis, s=5)
        axs[0, 0].set_xlabel('Muon θ (°)')
        axs[0, 0].set_ylabel('Muon φ (°)')
        axs[0, 0].set_title('Muon Angular Distribution')

        axs[0, 1].scatter(neutron_thetas, neutron_phis, s=5)
        axs[0, 1].set_xlabel('Neutron θ (°)')
        axs[0, 1].set_ylabel('Neutron φ (°)')
        axs[0, 1].set_title('Neutron Angular Distribution')

        axs[1, 0].scatter(rel_thetas, rel_phis, c=colors_rel, s=5)
        axs[1, 0].set_xlabel('Δθ (°)')
        axs[1, 0].set_ylabel('Δφ (°)')
        axs[1, 0].set_title('Relative Angular Difference')
        red_patch = mpatches.Patch(color='red', label='Other Event')
        blue_patch = mpatches.Patch(color='blue', label='Muon-Neutron Event')
        axs[1, 0].legend(handle)
        s=[red_patch, blue_patch], loc=('upper right')  # or adjust `loc` as needed

        bins = np.arange(0, max(separations) + 5, 5)
        axs[1, 1].hist(separations, bins=bins)
        axs[1, 1].set_xlabel('Separation Angle (°)')
        axs[1, 1].set_ylabel('Number of Events')
        axs[1, 1].set_title('Histogram of Separation Angles (5° bins)')

        plt.tight_layout()
        plt.savefig('../angle_distribution.pdf')


    def plot_length_inverse_hist_with_angle_curve(self):
        """
        Combine the first 10 events, collect all [length, energy, angle] rows,
        then:
          • Draw an energy histogram (bin = 0.1 MeV) with bar heights = counts
            on the left y-axis.
          • Compute the mean angle in each bin, and plot ⟨angle⟩ vs. bin_center
            on a secondary right y-axis as a line.
        Saves figure as: ./plots/energy_hist_angle_curve.png
        """

        # 1) Collect (energy, angle) from the first up to 10 events
        n_events = min(1000, self.event.nEntry)
        all_energy = []
        all_angle  = []
        all_length = []
        all_length_log = []
        all_length_inversesquare = []

        for idx in range(n_events):
            self.event.Jump(idx)
            data = self.event.edep_based_information()
            arr  = np.asarray(data)

            # If not N×3, skip
            if arr.ndim != 2 or arr.shape[1] != 3:
                print(
                    f"Warning: edep_based_information() for event {idx} "
                    f"returned shape {arr.shape}. Skipping."
                )
                continue

            # arr[:,1] = energy, arr[:,2] = angle
            all_length.extend(arr[:, 0].tolist())
            all_energy.extend(arr[:, 1].tolist())
            all_angle.extend(arr[:, 2].tolist())
            all_length_log.extend(np.log(arr[:, 0]).tolist())
            all_length_inversesquare.extend((arr[:,0]**(-2)).tolist())



        # 2) If no data was collected, bail out
        if len(all_energy) == 0:
            print("No valid (energy, angle) points found in the first 10 events; no histogram made.")
            return

        all_energy = np.asarray(all_energy)
        all_angle  = np.asarray(all_angle)
        all_length = np.asarray(all_length)
        all_length_inversesquare = np.asarray(all_length_inversesquare)
        all_length_log = np.asarray(all_length_log)
        #
        # # 3) Build bins of width 0.1 MeV
        # emin = all_length.min()
        # emax = all_energy.max()
        # # emin=0.5
        # emax=1.5

        #inverse case
        emin = 0
        emax =0.075


        bin_start = emin
        bin_stop = emax
        # bin_start = np.floor(emin * 10) / 10
        # bin_stop  = np.ceil(emax * 10) / 10
        bins = np.arange(bin_start, bin_stop , 0.005)
        counts, bin_edges = np.histogram(all_length_inversesquare, bins=bins)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        # 4) For each bin, compute angle‐mean
        means = []
        for i in range(len(bin_centers)):
            low  = bin_edges[i]
            high = bin_edges[i + 1]
            mask = (all_length_inversesquare >= low) & (all_length_inversesquare < high)

            if np.any(mask):
                angs = all_angle[mask]
                means.append(np.mean(angs))
            else:
                means.append(np.nan)

        # 5) Plot histogram on left y-axis
        fig, ax = plt.subplots(figsize=(8, 5), dpi=100)
        bars = ax.bar(
            bin_centers,
            counts,
            width=0.01,
            align="center",
            edgecolor="black",
            label="Counts"
        )
        ax.set_xlabel("inverse square length")
        ax.set_ylabel("number of points")
        ax.set_title(f"inverse square length histogram & ⟨angle⟩ curve (first {n_events} events)")

        # 6) Create a secondary y-axis for ⟨angle⟩
        ax2 = ax.twinx()
        ax2.plot(
            bin_centers,
            means,
            color="red",
            marker="o",
            linestyle="-",
            label="⟨angle⟩"
        )
        ax2.set_ylabel("mean angle")
        ax2.grid(False)  # turn off grid on the second axis

        # 7) Put a legend for the angle curve on the right
        ax2.legend(loc="upper right")

        # 8) Save to ./plots
        outdir = self.event.plotpath
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, "inversesquarelength_hist_angle_curve.png")
        fig.tight_layout()
        fig.savefig(outpath)
        plt.close(fig)
        print(f"  → saved combined histogram + ⟨angle⟩ curve to {outpath}")

    def plot_energy_vs_angle(self):

        # How many events to combine:
        n_events = min(10, self.event.nEntry)

        # Lists to accumulate every event’s points
        all_length = []
        all_energy = []
        all_angle  = []

        for idx in range(n_events):
            # 1) Jump into the idx-th event
            self.event.Jump(idx)

            # 2) Retrieve the raw data array (shape should be (N,3))
            data = self.event.edep_based_information()
            arr  = np.asarray(data)

            # 3) Verify it’s a 2D array with 3 columns; skip if not
            if arr.ndim != 2 or arr.shape[1] != 3:
                print(
                    f"Warning: edep_based_information() returned shape {arr.shape} "
                    f"for event {idx}. Skipping."
                )
                continue

            # 4) Unpack length/energy/angle columns
            length = arr[:, 0]
            energy = arr[:, 1]
            angle  = arr[:, 2]

            # 5) Append into our “big” lists
            all_length.extend(length.tolist())
            all_energy.extend(energy.tolist())
            all_angle.extend(angle.tolist())

        # If nothing was collected (all events skipped), warn & return
        if len(all_length) == 0:
            print("No valid points found in the first 10 events; no plot made.")
            return

        # Convert to numpy arrays for plotting
        all_length = np.asarray(all_length)
        all_energy = np.asarray(all_energy)
        all_angle  = np.asarray(all_angle)

        # 6) Single scatter: angle on X, energy on Y, color by length
        plt.figure(figsize=(7, 5), dpi=100)
        sc = plt.scatter(
            all_angle,
            all_energy,
            c=all_length,
            cmap="viridis",
            s=20,
            edgecolors="",
        )
        cbar = plt.colorbar(sc)
        cbar.set_label("length")
        plt.xlabel("angle")
        plt.ylabel("energy")
        plt.title(f"Combined energy vs. angle (first {n_events} events)")

        # 7) Ensure output directory exists and save
        outdir = self.event.plotpath
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, "energy_vs_angle_combined.png")
        plt.tight_layout()
        plt.savefig(outpath)
        plt.close()
        print(f"  → saved combined plot to {outpath}")


    def plot_length_vs_angle(self):

        # How many events to combine:
        n_events = min(10, self.event.nEntry)

        # Lists to accumulate every event’s points
        all_length = []
        all_energy = []
        all_angle  = []
        all_length_log = []
        all_length_squareinverse = []

        for idx in range(n_events):
            # 1) Jump into the idx-th event
            self.event.Jump(idx)

            # 2) Retrieve the raw data array (shape should be (N,3))
            data = self.event.edep_based_information()
            arr  = np.asarray(data)

            # 3) Verify it’s a 2D array with 3 columns; skip if not
            if arr.ndim != 2 or arr.shape[1] != 3:
                print(
                    f"Warning: edep_based_information() returned shape {arr.shape} "
                    f"for event {idx}. Skipping."
                )
                continue

            # 4) Unpack length/energy/angle columns
            length = arr[:, 0]
            energy = arr[:, 1]
            angle  = arr[:, 2]

            # 5) Append into our “big” lists
            all_length.extend(length.tolist())
            all_energy.extend(energy.tolist())
            all_angle.extend(angle.tolist())
            all_length_log.extend(log(length).tolist())
            all_length_squareinverse.extend((length**(-2).tolist))

        # If nothing was collected (all events skipped), warn & return
        if len(all_length) == 0:
            print("No valid points found in the first 10 events; no plot made.")
            return

        # Convert to numpy arrays for plotting
        all_length = np.asarray(all_length)
        all_energy = np.asarray(all_energy)
        all_angle  = np.asarray(all_angle)

        # 6) Single scatter: angle on X, energy on Y, color by length
        plt.figure(figsize=(7, 5), dpi=100)
        sc = plt.scatter(
            all_angle,
            all_energy,
            c=all_length,
            cmap="viridis",
            s=20,
            edgecolors="",
        )
        cbar = plt.colorbar(sc)
        cbar.set_label("length")
        plt.xlabel("angle")
        plt.ylabel("energy")
        plt.title(f"Combined energy vs. angle (first {n_events} events)")

        # 7) Ensure output directory exists and save
        outdir = self.event.plotpath
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, "energy_vs_angle_combined.png")
        plt.tight_layout()
        plt.savefig(outpath)
        plt.close()
        print(f"  → saved combined plot to {outpath}")





    #---------------------------------------------
    def DrawROOT(self, dim2d='yz', markerSize=0.2):
        import ROOT
        from ROOT import TH2F, TMarker, TCanvas, TLatex
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetMarkerStyle(24)
        ROOT.gStyle.SetMarkerSize(markerSize)
        colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kMagenta]
        nColor = len(colors)

        c1 = TCanvas("c1", self.event.vertex.GetReaction(), 800, 800)

        # canvas of 5m x 5m
        dummy = TH2F("dummy", "", 100, -5, 5, 100, -5, 5)
        dummy.GetXaxis().SetTitle('[m]')
        dummy.GetYaxis().SetTitle('[m]')
        dummy.Draw()

        mm2m = 0.001
        txts = []
        markers = []
        txtX = 0.15
        txtY = 0.85
        txtSize = 0.03
        # nDepo = 0
        for i, track in enumerate(self.event.tracks):
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            color = colors[ancestor % nColor]


            if ancestor == i:
                txt = TLatex()
                txt = txt.DrawLatexNDC(txtX, txtY, track.GetName())
                txt.SetTextColor(color)
                txt.SetTextSize(txtSize)
                txts.append(txt)
                txtY -= txtSize

            for j, di in enumerate(depoList):
                depo = self.event.depos[di]
                x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                mapping = {'x': x, 'y': y, 'z': z}

                m = TMarker(mapping[dim2d[1]], mapping[dim2d[0]], 24)
                m.SetMarkerColor(color)

                m.Draw()
                markers.append(m)
                # nDepo += 1

        # print('depo points drawn: ', nDepo, '| total depo: ', self.event.depos.size)

        ROOT.gPad.Update()
        return c1



    def effeciency_plot(self, markerSize=3, cmap='viridis'):
        # --- prepare data & figure with shared axes (item 4) ---
        self.event.Jump(self.event.currentEntry)
        self.Collect(0)

        fig, axes = plt.subplots(
            nrows=2, ncols=3,
            figsize=(9*3, 9*2), dpi=500,
            sharex=True, sharey=True
        )
        ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()

        # --- 2) ax1: Y vs Z colored by time + colorbar (item 2) ---
        sc = ax1.scatter(
            self.zz, self.yy,
            c=self.tt,
            cmap=cmap,
            s=15,
            norm=matplotlib.colors.LogNorm(),
            edgecolors='none',
            linewidths=0,
            zorder=1
        )
        ax1.set_xlabel('Z (m)')
        ax1.set_ylabel('Y (m)')
        ax1.set_title('Y vs Z, colored by time')
        fig.colorbar(sc, ax=ax1, label='Time (ns)')

        # capture the common limits
        x0, x1 = ax1.get_xlim()
        y0, y1 = ax1.get_ylim()

        # --- 1) ax4: fill the previously empty slot with a PDG‐ID histogram (item 1) ---
        pdg_codes = [
            p.GetPDGCode()
            for p in self.event.vertex.Particles
            if p.GetTrackId() >= 0
        ]
        ax4.hist(pdg_codes, bins=np.arange(min(pdg_codes), max(pdg_codes) + 2) - 0.5)
        ax4.set_title("PDG ID distribution")
        ax4.set_xlabel("PDG Code")
        ax4.set_ylabel("Count")

        # --- compute and mask your cone data as before ---
        real_neutron = np.array(self.event.reconstructing_direction())
        filtered_neutron_1 = np.array(self.event.coneexclusion(10, 0))
        filtered_neutron_2 = np.array(self.event.coneexclusion(20, 0))
        print("the real_neutron is :",real_neutron)
        print("the filtered_neutron_1 is :",filtered_neutron_1)
        print("the filtered_neutron_2 is :",filtered_neutron_2)

        if len(real_neutron) == 0 or len(filtered_neutron_1) == 0:
            print("empty data, skipping plots")
            return

        y_r, z_r = real_neutron[:, 2], real_neutron[:, 3]
        y_f1, z_f1 = filtered_neutron_1[:, 2], filtered_neutron_1[:, 3]
        y_f2, z_f2 = filtered_neutron_2[:, 2], filtered_neutron_2[:, 3]

        # --- 3) ax3 & ax5 & ax6: the 3 neutron‐deposition scatter plots ---
        ax3.scatter(z_r, y_r, s=15, edgecolors='none')
        ax3.set_title("Real neutron deposition")
        ax3.set_ylabel("Y(m)")
        ax3.set_xlabel("Z(m)")

        ax5.scatter(z_f1, y_f1, s=15, edgecolors='none')
        ax5.set_title("Filtered (10°)")
        ax5.set_xlabel('Z (m)')
        ax5.set_ylabel('Y (m)')

        ax6.scatter(z_f2, y_f2, s=15, edgecolors='none')
        ax6.set_title("Filtered (20°)")

        ax2.set_xlabel('Z (m)')
        ax2.set_ylabel('Y (m)')
        # --- 5) ax2: arrows + true legend via Patch handles (item 5) ---
        mask = self.ee >= 0.5
        ax2.scatter(
            self.zz[mask], self.yy[mask],
            c=self.cc[mask],
            s=markerSize,
            cmap='tab20'
        )
        # draw neutron & muon arrows
        comp = {'x': 0, 'y': 1, 'z': 2}
        dn = self.event.read_neutron_direction()
        dm = self.event.read_muon_direction()
        for vec, color, label in ((dn, 'blue', 'neutron'), (dm, 'red', 'muon')):
            v = np.array([vec[comp['y']], vec[comp['z']]])
            v = v / np.linalg.norm(v)
            ax2.arrow(0, 0, v[0], v[1],
                      color=color, length_includes_head=True,
                      head_width=0.1, head_length=0.1)
        # legend via patches
        handles = [
            mpatches.Patch(color='blue', label='neutron'),
            mpatches.Patch(color='red', label='muon'),
        ]
        ax2.legend(handles=handles, loc='upper right', fontsize=8)
        ax2.set_title("Ancestor IDs & directions")
        zdata = self.zz[mask]
        ydata = self.yy[mask]
        # add a 10% padding, for instance
        zpad = 0.1 * (zdata.max() - zdata.min())
        ypad = 0.1 * (ydata.max() - ydata.min())

        ax2.set_xlim(zdata.min() - zpad, zdata.max() + zpad)
        ax2.set_ylim(ydata.min() - ypad, ydata.max() + ypad)
        ax2.set_xlabel('Z (m)')
        ax2.set_ylabel('Y (m)')
        # --- 6) remove plt.clf() (item 6) ---
        # re-apply the shared limits to all panels once:
        for ax in (ax3, ax4, ax5, ax6):
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
        for ax in axes.flatten():
            ax.tick_params(labelbottom=True, labelleft=True)
        fig.savefig(f"{self.event.plotpath}{self.event.currentEntry}_cone_exclusion.pdf")
        plt.close(fig)



if __name__ == "__main__":
    event = Event(sys.argv[1])
    p = Plotter(event)
    #p.Next()
    #p.Draw('yz')
    # c1 = p.DrawROOT('xz', 0.2)
    # input('press a key to continue ...')
    for entry in range(20):
        p.event.currentEntry = entry
    #    p.Draw('yz')
        p.effeciency_plot()
    # p.plot_length_inverse_hist_with_angle_curve()
