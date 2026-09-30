# Template matching by normalised cross-correlation. The template is a matrix of samples.
SeedRandom[3]; scene = RandomImage[1, {40, 30}];
tmpl = ImageData[scene][[11 ;; 17, 21 ;; 27]];
ncc = ImageData[ImageCorrelate[scene, tmpl, "NormalizedCrossCorrelation"]];
{Position[ncc, Max[ncc]], Max[ncc]}
data = ImageData[scene]; data[[22 ;; 30, 1 ;; 9]] = 1.;
bright = Image[0.4 data + 0.3];
raw = ImageData[ImageCorrelate[bright, tmpl]];
Position[raw, Max[raw]][[1]]
ncc2 = ImageData[ImageCorrelate[bright, tmpl, "NormalizedCrossCorrelation"]];
Position[ncc2, Max[ncc2]]
