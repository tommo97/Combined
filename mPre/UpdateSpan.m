function Span = UpdateSpan(Span)
Span.y = [];
Span = SetBladeData(Span);
Span.DistPanel = PanelDistButtonsParam(Span.DistPanel);
Span.y = Span.DistPanel.x;


Span = SetBladeData(Span);



