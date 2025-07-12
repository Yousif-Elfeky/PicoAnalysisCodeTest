#ifndef CUTSCONFIG_H
#define CUTSCONFIG_H

// Central location for all analysis cut values.
// Adjust for different energies by editing these constants or
// generating them from an external config later.

#include "StPicoEvent/StPicoEvent.h"

namespace HFCuts {
    // -------- Event cuts (common) --------
    namespace Event {
        constexpr float VzMax      = 60.f;   // |Vz| (cm) for collider; ignore for FXT with special window below
        constexpr float VzDiffMax  = 6.f;    // |Vz - VzVPD|
        // Fixed-target window along beamline (cm)
        constexpr float VzMinFXT   = 185.f;
        constexpr float VzMaxFXT   = 235.f;
    }

    // -------- Track quality --------
    namespace Track {
        constexpr int   nHitsFitMin    = 15;   // FXT: 15
        constexpr float nHitsFracMin   = 0.52;
        constexpr float dcaMax         = 3.0;  // FXT: 3.0
        constexpr float ptMin          = 0.1f; // FXT: 0.1
    }

    // -------- PID windows --------
    namespace PID {
        constexpr float betaDiffPiMax = 0.03; // |1/β - 1/β_expected| cut for pions
        constexpr float betaDiffKMax  = 0.03; // for kaons
        // TPC nσ
        constexpr float nSigmaE   = 3.0;
        constexpr float nSigmaPi  = 2.0;
        constexpr float nSigmaK   = 2.0;
        // TOF 1/beta cut (absolute)
        constexpr float oneOverBetaMax = 0.05; // 19 GeV, loosen to 0.04 at FXT
        // E/p window for electrons (BEMC)
        constexpr float eopMin = 0.8;
        constexpr float eopMax = 1.2;
    }

    // -------- helper to evaluate event cuts --------
    inline bool passEventCuts(const StPicoEvent* evt) {
        if (!evt) return false;
        float vz = evt->primaryVertex().z();
        float vzVpd = evt->vzVpd();
        return std::fabs(vz) < Event::VzMax && std::fabs(vz - vzVpd) < Event::VzDiffMax;
    }

    // -------- Physics observable useful ranges --------
    namespace Obs {
        constexpr float JPsiMassMin = 2.0;
        constexpr float JPsiMassMax = 4.0;
        constexpr float D0MassMin   = 1.6;
        constexpr float D0MassMax   = 2.1;
    }
}

#endif
