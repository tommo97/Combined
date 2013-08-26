function Span = UpdateSpan(Span)
Span.y = [];
Span = SetBladeData(Span);
Span.DistPanel = HalfPanelDistButtonsParam(Span.DistPanel);
Span.y = Span.DistPanel.x;


Span = SetBladeData(Span);



