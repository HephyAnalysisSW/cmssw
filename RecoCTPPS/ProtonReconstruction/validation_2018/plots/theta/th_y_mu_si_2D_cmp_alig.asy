import root;
import pad_layout;

string topDir = "../../data_eos/";

string stream = "ZeroBias";

string dataset = "fill_7334/xangle_140_beta_0.30";

string alignments[], a_labels[];
alignments.push("2018_11_02.3"); a_labels.push("2018-11-02.3");

string cols[], c_labels[];
cols.push("si_rp3_mu_arm0"); c_labels.push("multi(L) vs.~single(RP3)");
cols.push("si_rp23_mu_arm0"); c_labels.push("multi(L) vs.~single(RP23)");
cols.push("si_rp103_mu_arm1"); c_labels.push("multi(R) vs.~single(RP103)");
cols.push("si_rp123_mu_arm1"); c_labels.push("multi(R) vs.~single(RP123)");

TH2_palette = Gradient(blue, heavygreen, yellow, red);

//----------------------------------------------------------------------------------------------------

NewPad(false);
label("\vbox{\hbox{stream: " + stream + "}\hbox{dataset: " + replace(dataset, "_", "\_") + "}}");

NewRow();

NewPad(false);
for (int ci : cols.keys)
	NewPadLabel(c_labels[ci]);

for (int ai : alignments.keys)
{
	NewRow();

	NewPadLabel(a_labels[ai]);

	for (int ci : cols.keys)
	{
		NewPad("$\th^*_y (\rm single)\ung{\mu rad}$", "$\th^*_y (\rm multi)\ung{\mu rad}$", axesAbove=true);
		scale(Linear, Linear, Log);

		string f = topDir + dataset + "/" + stream + "/alignment_" + alignments[ai] + "/output.root";
		string on = "singleMultiCorrelationPlots/" + cols[ci] + "/h2_th_y_mu_vs_th_y_si";
		
		RootObject obj = RootGetObject(f, on, error=true);
		if (!obj.valid)
			continue;
		
		draw(scale(1e6, 1e6), obj);

		draw((-200, -200)--(+200, +200), black+1pt);

		limits((-200, -200), (+200, +200), Crop);

		xaxis(YEquals(0, false), black+1pt, above=true);
		yaxis(XEquals(0, false), black+1pt, above=true);
	}
}

//----------------------------------------------------------------------------------------------------

GShipout(hSkip=1mm, vSkip=0mm);
