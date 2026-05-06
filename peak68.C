// peak68.C
// Compute 68.27% containment interval around a user-input peak
// Works for 1D histograms (TH1)

#include <iostream>
#include <vector>
#include <cmath>

double TARGET_FRAC = 0.682689492137;  // 1 sigma

namespace {

double integrate_hist(TH1* h, double x_min, double x_max, int nsteps = 2000) {
    if (x_max <= x_min)
        return 0.0;
    double step = (x_max - x_min) / nsteps;
    double sum = 0.0;
    for (int i = 0; i <= nsteps; ++i) {
        double x = x_min + i * step;
        double val = h->Interpolate(x);
        if (val < 0.0)
            val = 0.0;
        double weight = (i == 0 || i == nsteps) ? 0.5 : 1.0;
        sum += weight * val;
    }
    return sum * step;
}

} // anonymous namespace

void compute_interval_TH1(TH1* h, double x_peak) {

    int nb = h->GetNbinsX();
    if (nb <= 0) {
        std::cout << "Histogram has no bins.\n";
        return;
    }

    // cache bin centers and non-negative contents
    std::vector<double> centers(nb);
    std::vector<double> contents(nb);
    double total_area = 0.0;
    for (int i = 1; i <= nb; ++i) {
        centers[i-1] = h->GetBinCenter(i);
        double c = h->GetBinContent(i);
        if (c < 0) c = 0;
        contents[i-1] = c;
        total_area += c;
    }

    if (total_area <= 0) {
        std::cout << "Histogram empty.\n";
        return;
    }

    double target_area = TARGET_FRAC * total_area;

    int peak_bin = h->GetXaxis()->FindBin(x_peak) - 1;
    peak_bin = std::max(0, std::min(nb - 1, peak_bin));

    auto find_left_cross = [&](double level) {
        int i = peak_bin;
        while (i > 0 && contents[i] >= level)
            --i;

        if (contents[i] >= level)
            return h->GetXaxis()->GetXmin();

        int j = i + 1;
        double x1 = centers[i];
        double y1 = contents[i];
        double x2 = centers[j];
        double y2 = contents[j];
        if (std::fabs(y2 - y1) < 1e-9)
            return x1;
        double t = (level - y1) / (y2 - y1);
        return x1 + t * (x2 - x1);
    };

    auto find_right_cross = [&](double level) {
        int i = peak_bin;
        while (i < nb - 1 && contents[i] >= level)
            ++i;

        if (contents[i] >= level)
            return h->GetXaxis()->GetXmax();

        int j = i - 1;
        double x1 = centers[j];
        double y1 = contents[j];
        double x2 = centers[i];
        double y2 = contents[i];
        if (std::fabs(y2 - y1) < 1e-9)
            return x2;
        double t = (level - y1) / (y2 - y1);
        return x2 + t * (x1 - x2);
    };

    double max_height = 0.0;
    for (double c : contents)
        if (c > max_height)
            max_height = c;

    double lo = 0.0;
    double hi = max_height;

    double best_level = 0.0;
    double best_x_low = h->GetXaxis()->GetXmin();
    double best_x_high = h->GetXaxis()->GetXmax();

    for (int iter = 0; iter < 40; ++iter) {
        double mid = 0.5 * (lo + hi);
        double x_low = find_left_cross(mid);
        double x_high = find_right_cross(mid);
        if (x_high <= x_low) {
            hi = mid;
            continue;
        }

        double area = integrate_hist(h, x_low, x_high, 2000);

        if (area >= target_area) {
            best_level = mid;
            best_x_low = x_low;
            best_x_high = x_high;
            lo = mid;  // raise lower bound to tighten interval
        } else {
            hi = mid;
        }
    }

    double err_low  = x_peak - best_x_low;
    double err_high = best_x_high - x_peak;
    double interval_area = integrate_hist(h, best_x_low, best_x_high, 5000);

    std::cout << "\n=== Result for " << h->GetName() << " ===\n";
    std::cout << "Peak x = " << x_peak << "\n";
    std::cout << "Interval = [" << best_x_low << ", " << best_x_high << "]\n";
    std::cout << "Error bars: -" << err_low << "  +" << err_high << "\n";
    if (best_level > 0)
        std::cout << "Boundary height (common y) = " << best_level << "\n";
    std::cout << "Included fraction = " << interval_area / total_area << "\n";
}


void peak68() {

    // assume file already opened, or open here
    TFile *f = TFile::Open("comparison_with_lepton_neutron_event0.2.root");

    // get your canvas
    TCanvas *c =
        (TCanvas*)f->Get("first_reconstruction_result_on_neutrino");

    c->Draw();

    // loop over objects drawn in canvas
    auto prims = c->GetListOfPrimitives();

    for (int i=0;i<prims->GetSize();i++) {

        TObject* obj = prims->At(i);

        if (obj->InheritsFrom("TH1")) {

            TH1* h = (TH1*)obj;

            std::cout << "\nHistogram: " << h->GetName() << "\n";

            double x_peak;
            std::cout << "Enter peak x value: ";
            std::cin >> x_peak;

            compute_interval_TH1(h, x_peak);
        }
    }
}
