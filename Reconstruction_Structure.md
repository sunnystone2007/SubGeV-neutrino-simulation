From printtracks in event read energy deposition, pdg code, x,y,z coordinate  

add a "selectneutronevent" to make sure there is one neutron coming out of first vertex, if not,  select the one with smallest pdg code  

add a "neutron_direction" function in event, input the information in printtracks, print out the real nuetron direction, calculated neutron direction  
and the cos angle between two angles.    



    def selectneutronevent(self):
        trig = 0
        for i in range(self.tracks.size):
            track = self.tracks[i]

            pdg_original = track.GetPDGCode()
            pdg = track.GetPDGCode()  # Get the PDG code of the particle
            track_energy = track.energy.get('depoTotal', 0)
            # Trace back to the parent particle until we find a neutron (PDG code 2112)
            pdg, track = self.loopover(pdg, track)
            ParentId = track.GetParentId()
            while pdg == 2112 and ParentId != -1:
                track = self.tracks[ParentId]
                pdg = track.GetPDGCode()
                pdg, track = self.loopover(pdg, track)
                ParentId = track.GetParentId()
            if pdg == 2112 and ParentId == -1 and track_energy>0.5 and pdg_original!=2112:
                trig+=1
        if trig==1:
            print("this",self.currentEntry,"th event has one neutron neutrino interaction, it is a good event")
        elif trig!=0 and trig!=1:
            print("this",self.currentEntry,"th event has",trig,"neutron neutrino interaction, it is a BAD event")


      def PrintTracks(self, start=0, stop=-1):
        #print(f"{self.tracks.size} trajectories stored", )

        #print(f"{'pdg':>8}{'name':>8}{'trkId':>6}{'parId':>6}{'acId':>6}{'KE':>10}{'selfDepo':>10}{'allDepo':>10}")
        #print(f"{'':>8}{'':>8}{'':>6}{'':>6}{'':>6}{'[MeV]':>10}{'[MeV]':>10}{'[MeV]':>10}")
        #print('-'*(8+8+6+6+6+10+10+10))

        neutrontrkId = -1
        neutronKE = -1
        numbers=self.select_the_right_track()
        for track in self.tracks[numbers]:
            pdg = track.GetPDGCode()
            print("the particle is:",pdg)
            name = track.GetName()
            trkId = track.GetTrackId()
            parId = track.GetParentId()
            mom = track.GetInitialMomentum()
            # print("the momentum is :",mom)
            px = track.GetInitialMomentum().X()
            py = track.GetInitialMomentum().Y()
            pz = track.GetInitialMomentum().Z()
            pE = track.GetInitialMomentum().E()
            pM = track.GetInitialMomentum().M()
            for point in track.Points:
               x = point.GetPosition().X()
               y = point.GetPosition().Y()
               z = point.GetPosition().Z()
               t = point.GetPosition().T()
            print("the x,y,z,t coordinate of the track:",x,y,z,t)
            #print("the 3 momentum are:",px,py,pz,"the energy is:",pE,"the mass is:",pM,)
            mass = mom.M()
            KE = mom.E() - mass
            print("the kenitic energy is: ",KE)
            p_square=(px**2+py**2+pz**2)**0.5
            direction_vector=[]
            direction_vector.append(px/p_square)
            direction_vector.append(py/p_square)
            direction_vector.append(pz/p_square)
            if KE>=2:
                print(direction_vector)
            else:
                print("direction undetermined")
            ancestor = track.association['ancestor']
            selfDepo = track.energy['depoTotal']
            allDepo = self.GetEnergyDepoWithDesendents(trkId)
            #print(f"{pdg:>8d}{name:>8s}{trkId:>6d}{parId:>6d}{ancestor:>6d}{KE:>10.2f}{selfDepo:>10.2f}{allDepo:>10.2f}")
            # Take the last neutron and check its KE and capture time
            if pdg == 2112:
                neutrontrkId = trkId
                neutronKE = KE

        #print('-'*(8+8+6+6+6+10+10+10))
        return [neutrontrkId, neutronKE]
