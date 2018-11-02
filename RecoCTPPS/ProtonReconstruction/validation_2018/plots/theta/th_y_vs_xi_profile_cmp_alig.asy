import root;
import pad_layout;

string topDir = "../../data_eos/";

string dataset = "fill_7333/xangle_140_beta_0.30";

string stream = "ZeroBias";

string alignments[], a_labels[];
//alignments.push("2018_11_02.1"); a_labels.push("2018-11-02.1");
//alignments.push("2018_11_02.2"); a_labels.push("2018-11-02.2");
alignments.push("2018_11_02.3"); a_labels.push("2018-11-02.3");

string cols[], c_labels[];
cols.push("arm0"); c_labels.push("sector 45 (L)");
cols.push("arm1"); c_labels.push("sector 56 (R)");

TH2_palette = Gradient(blue, heavygreen, yellow, red);

//----------------------------------------------------------------------------------------------------

NewPad(false);
AddToLegend("stream = " + stream);
AddToLegend("dataset = " + replace(dataset, "_", "\_"));
AttachLegend();

for (int ci : cols.keys)
	NewPadLabel(c_labels[ci]);

for (int ai : alignments.keys)
{
	NewRow();

	NewPadLabel(a_labels[ai]);

	for (int ci : cols.keys)
	{
		NewPad("$\xi_{\rm multi}$", "mean of $\th^*_y\ung{\mu rad}$");

		string f = topDir + dataset + "/" + stream + "/alignment_" + alignments[ai] + "/do_fits.root";
		string on = "multiRPPlots/" + cols[ci] + "/p_th_y_vs_xi";

		RootObject hist = RootGetObject(f, on);
		RootObject fit = RootGetObject(f, on + "|ff");

		draw(scale(1., 1e6), hist, "d0,eb", blue);
		draw(scale(1., 1e6), fit, "def", red+1pt);

		limits((0.00, -300), (0.20, +300), Crop);
	}
}

GShipout(hSkip=0mm, vSkip=0mm);
