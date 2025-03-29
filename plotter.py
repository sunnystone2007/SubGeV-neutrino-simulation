from event import Event
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

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

        self.USER_COLORS = ['black', 'red', 'blue', 'magenta']

        nColor = len(self.USER_COLORS)
        mm2m = 0.001
        mm2cm = 0.1
        for i, track in enumerate(self.event.tracks):
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            for di in depoList:
                depo = self.event.depos[di]
                x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                t = (depo.GetStart().T() + depo.GetStop().T()) / 2 # ns
                e = depo.GetEnergyDeposit() # MeV
                l = depo.GetTrackLength() *mm2cm # most are 0.5 cm
                if self.event.ChargeBirksLaw(e, l) > dQthreshold: # detector threshold
                    self.xx = np.append(self.xx, x)
                    self.yy = np.append(self.yy, y)
                    self.zz = np.append(self.zz, z)
                    self.tt = np.append(self.tt, t)
                    self.ee = np.append(self.ee, e)
                    self.ll = np.append(self.ll, l)
                    self.cc = np.append(self.cc, self.USER_COLORS[ancestor % nColor])

    def Jump(self, entryNo, dQthreshold):
        self.event.Jump(entryNo)
        self.Collect(dQthreshold)

    def Next(self, dQthreshold):
        self.event.Next()
        self.Collect(dQthreshold)

    def Prev(self, dQthreshold):
        self.event.Prev()
        self.Collect(dQthreshold)

    def Draw(self, axis='yz', value='time', energy= 'GeV', markerSize=3, cmap='jet', vmax=2000):
        # particle, timing, dE/dx
        mapping = {'x': self.xx, 'y': self.yy, 'z': self.zz}

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5*2, 4), dpi=100)
        fig.suptitle(self.event.vertex.GetReaction())
        cb_ax = fig.add_axes([.94,.124,.02,.754])

        # particle plot
        ax1.scatter(mapping[axis[1]], mapping[axis[0]], c=self.cc, s=markerSize)
        if value == 'time':
            # timing plot
            plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.tt, cmap=cmap, vmin=1, vmax=vmax, norm=matplotlib.colors.LogNorm(), s=markerSize)
            cb_ax.set_xlabel('ns')

        elif value == 'charge':
            # charge plot
            if energy == 'GeV':
                plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.ee*2, cmap=cmap, vmax=4, s=markerSize)
            elif energy == 'MeV':
                #plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=np.divide(self.ee, self.ll), cmap=cmap, vmax=4, s=markerSize)
                plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.ee*2, cmap=cmap, vmax=4, s=markerSize)
            cb_ax.set_xlabel('MeV/cm')


        fig.colorbar(plot_12, orientation='vertical', cax=cb_ax)
        # fig.tight_layout()

        for ax in (ax1, ax2):
            if self.event.evgen == 'Genie':
                ax.set_ylim(-4, 4)
                ax.set_xlim(-2, 8)
                xpos = -1.8
                ypos = 3.6
                interval = 0.35
            elif self.event.evgen == 'Marley':
                # Normally Marley use this scale
                #ax.set_ylim(-1, 1.5)
                #ax.set_xlim(-1, 1.5)
                #interval = 0.1
                # Debug neutron depo
                #ax.set_ylim(-200, 200)
                #ax.set_xlim(-200, 200)
                #interval = 16
                #xpos = -195
                #ypos = 190
                # n-capture scale
                #xpos = -0.45
                #ypos = 0.4
                #ax.set_ylim(-0.5, 0.5)
                #ax.set_xlim(-0.5, 0.5)
                xpos = -0.95
                ypos = 0.9
                ax.set_ylim(-1, 1)
                ax.set_xlim(-1, 0.5)
                interval = 0.1
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
            color = self.USER_COLORS[(i-countnegId) % nColor]
            # pdg = particle.GetPDGCode()
            # name = particle.GetName()
            # trkId = particle.GetTrackId()
            mom = particle.GetMomentum()
            KE = mom.E() - mom.M()
            name = '%s: %.1f MeV' % (name, KE)

            ax1.text(xpos, ypos, name, color=color)
            ypos -= interval

        #ax2.plot()
        fig.savefig(self.event.plotpath + '/particle_%s_%s_%s_evt_%d.pdf' % (value, axis, energy, self.event.currentEntry) )
        plt.clf() # important to clear figure
        plt.close()
        #plt.show()

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

if __name__ == "__main__":
    event = Event(sys.argv[1])
    p = Plotter(event)
    #p.Next()
    #p.Draw('yz')
    # c1 = p.DrawROOT('xz', 0.2)
    # input('press a key to continue ...')