function Chord = UpdateFoil(Chord)
Chord.DistPanel = HalfPanelDistButtonsParam(Chord.DistPanel);
Chord.DistPanel.x = fliplr(1-Chord.DistPanel.x);
Chord.Tip = SetFoilData(Chord.Tip,Chord.DistPanel);
Chord.Root = SetFoilData(Chord.Root,Chord.DistPanel);



