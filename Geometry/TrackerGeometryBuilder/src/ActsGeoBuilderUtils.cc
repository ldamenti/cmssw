#include "Geometry/TrackerGeometryBuilder/interface/ActsGeoBuilderUtils.h"

  std::vector<std::shared_ptr<Acts::Surface>> ActsGeoBuilderUtils::SelectActiveSurfaces_PhaseII(const KdtSurfacesDim2Bin100& surfaces, const std::string Layer_name){

    Acts::RangeXD<2, double, std::array> Range({0, 0}, {0, 0});

    // TBPX
    if(Layer_name=="TBPX_L1"){
      Range = Acts::RangeXD<2, double, std::array>({-200, 30}, {200, 43});
    }
    if(Layer_name=="TBPX_L2"){
      Range = Acts::RangeXD<2, double, std::array>({-200, 57}, {200, 73});
    }
    if(Layer_name=="TBPX_L3"){
      Range = Acts::RangeXD<2, double, std::array>({-200, 100}, {200, 117});
    }
    if(Layer_name=="TBPX_L4"){
      Range = Acts::RangeXD<2, double, std::array>({-200, 142}, {200, 162});
    }
    // TFPX positive
    if(Layer_name=="TFPX_D1_pos"){
      Range = Acts::RangeXD<2, double, std::array>({238, 30}, {266, 165});
    }
    if(Layer_name=="TFPX_D2_pos"){
      Range = Acts::RangeXD<2, double, std::array>({308, 30}, {336, 165});
    }
    if(Layer_name=="TFPX_D3_pos"){
      Range = Acts::RangeXD<2, double, std::array>({400, 30}, {420, 165});
    }
    if(Layer_name=="TFPX_D4_pos"){
      Range = Acts::RangeXD<2, double, std::array>({511, 30}, {539, 165});
    }
    if(Layer_name=="TFPX_D5_pos"){
      Range = Acts::RangeXD<2, double, std::array>({658, 30}, {679, 165});
    }
    if(Layer_name=="TFPX_D6_pos"){
      Range = Acts::RangeXD<2, double, std::array>({826, 30}, {854, 165});
    }
    if(Layer_name=="TFPX_D7_pos"){
      Range = Acts::RangeXD<2, double, std::array>({1099, 40}, {1120, 165});
    }
    if(Layer_name=="TFPX_D8_pos"){
      Range = Acts::RangeXD<2, double, std::array>({1386, 50}, {1407, 165});
    }
    // TFPX negative
    if(Layer_name=="TFPX_D1_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-266, 30}, {-238, 165});
    }
    if(Layer_name=="TFPX_D2_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-336, 30}, {-308, 165});
    }
    if(Layer_name=="TFPX_D3_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-420, 30}, {-400, 165});
    }
    if(Layer_name=="TFPX_D4_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-539, 30}, {-511, 165});
    }
    if(Layer_name=="TFPX_D5_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-679, 30}, {-658, 165});
    }
    if(Layer_name=="TFPX_D6_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-854, 30}, {-826, 165});
    }
    if(Layer_name=="TFPX_D7_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-1120, 40}, {-1099, 165});
    }
    if(Layer_name=="TFPX_D8_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-1407, 50}, {-1386, 165});
    }
    // TEPX positive
    if(Layer_name=="TEPX_D1_pos"){
      Range = Acts::RangeXD<2, double, std::array>({1743, 63}, {1757, 260});
    }
    if(Layer_name=="TEPX_D2_pos"){
      Range = Acts::RangeXD<2, double, std::array>({2000, 73}, {2016, 260});
    }
    if(Layer_name=="TEPX_D3_pos"){
      Range = Acts::RangeXD<2, double, std::array>({2290, 83}, {2320, 260});
    }
    if(Layer_name=="TEPX_D4_pos"){
      Range = Acts::RangeXD<2, double, std::array>({2630, 93}, {2660, 260});
    }
    // TEPX negative
    if(Layer_name=="TEPX_D1_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-1757, 63}, {-1743, 260});
    }
    if(Layer_name=="TEPX_D2_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-2016, 73}, {-2000, 260});
    }
    if(Layer_name=="TEPX_D3_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-2320, 83}, {-2290, 260});
    }
    if(Layer_name=="TEPX_D4_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-2660, 93}, {-2630, 260});
    }
    // TEDD positive
    if(Layer_name=="TEDD_D1_pos"){
      Range = Acts::RangeXD<2, double, std::array>({1281, 222}, {1340, 1100});
    }
    if(Layer_name=="TEDD_D2_pos"){
      Range = Acts::RangeXD<2, double, std::array>({1500, 222}, {1580, 1100});
    }
    if(Layer_name=="TEDD_D3_pos"){
      Range = Acts::RangeXD<2, double, std::array>({1815, 300}, {1890, 1100});
    }
    if(Layer_name=="TEDD_D4_pos"){
      Range = Acts::RangeXD<2, double, std::array>({2160, 300}, {2250, 1100});
    }
    if(Layer_name=="TEDD_D5_pos"){
      Range = Acts::RangeXD<2, double, std::array>({2600, 300}, {2690, 1100});
    }
    // TEDD negative
    if(Layer_name=="TEDD_D1_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-1340, 222}, {-1281, 1100});
    }
    if(Layer_name=="TEDD_D2_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-1580, 222}, {-1500, 1100});
    }
    if(Layer_name=="TEDD_D3_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-1890, 300}, {-1815, 1100});
    }
    if(Layer_name=="TEDD_D4_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-2250, 300}, {-2160, 1100});
    }
    if(Layer_name=="TEDD_D5_neg"){
      Range = Acts::RangeXD<2, double, std::array>({-2690, 300}, {-2600, 1100});
    }
    // TB2S
    if(Layer_name=="TB2S_L1"){
      Range = Acts::RangeXD<2, double, std::array>({-1230, 660}, {1230, 680});
    }
    if(Layer_name=="TB2S_L2"){
      Range = Acts::RangeXD<2, double, std::array>({-1230, 690}, {1230, 715});
    }
    if(Layer_name=="TB2S_L3"){
      Range = Acts::RangeXD<2, double, std::array>({-1230, 835}, {1230, 850});
    }
    if(Layer_name=="TB2S_L4"){
      Range = Acts::RangeXD<2, double, std::array>({-1230, 870}, {1230, 890});
    }
    if(Layer_name=="TB2S_L5"){
      Range = Acts::RangeXD<2, double, std::array>({-1230, 1060}, {1230, 1075});
    }
    if(Layer_name=="TB2S_L6"){
      Range = Acts::RangeXD<2, double, std::array>({-1230, 1090}, {1230, 1115});
    }
    // TBPS barrel
    if(Layer_name=="TBPS_L1"){
      Range = Acts::RangeXD<2, double, std::array>({-160, 210}, {160, 236});
    }
    if(Layer_name=="TBPS_L2"){
      Range = Acts::RangeXD<2, double, std::array>({-160, 238}, {160, 260});
    }
    if(Layer_name=="TBPS_L3"){
      Range = Acts::RangeXD<2, double, std::array>({-259, 340}, {259, 360});
    }
    if(Layer_name=="TBPS_L4"){
      Range = Acts::RangeXD<2, double, std::array>({-259, 367}, {259, 383});
    }
    if(Layer_name=="TBPS_L5"){
      Range = Acts::RangeXD<2, double, std::array>({-350, 490}, {350, 512});
    }
    if(Layer_name=="TBPS_L6"){
      Range = Acts::RangeXD<2, double, std::array>({-350, 520}, {350, 535});
    }
    // TBPS tilted L1 pos
    if(Layer_name=="TBPS_PosTilted_L1D1"){
      Range = Acts::RangeXD<2, double, std::array>({147, 213}, {196, 300});
    }
    if(Layer_name=="TBPS_PosTilted_L1D2"){
      Range = Acts::RangeXD<2, double, std::array>({197, 213}, {245, 300});
    }
    if(Layer_name=="TBPS_PosTilted_L1D3"){
      Range = Acts::RangeXD<2, double, std::array>({246, 213}, {294, 300});
    }
    if(Layer_name=="TBPS_PosTilted_L1D4"){
      Range = Acts::RangeXD<2, double, std::array>({295, 213}, {343, 300});
    }
    if(Layer_name=="TBPS_PosTilted_L1D5"){
      Range = Acts::RangeXD<2, double, std::array>({350, 213}, {406, 300});
    }
    if(Layer_name=="TBPS_PosTilted_L1D6"){
      Range = Acts::RangeXD<2, double, std::array>({420, 213}, {476, 300});
    }
    if(Layer_name=="TBPS_PosTilted_L1D7"){
      Range = Acts::RangeXD<2, double, std::array>({504, 213}, {553, 300});
    }
    if(Layer_name=="TBPS_PosTilted_L1D8"){
      Range = Acts::RangeXD<2, double, std::array>({595, 213}, {630, 300});
    }
    if(Layer_name=="TBPS_PosTilted_L1D9"){
      Range = Acts::RangeXD<2, double, std::array>({700, 213}, {742, 300});
    }
    if(Layer_name=="TBPS_PosTilted_L1D10"){
      Range = Acts::RangeXD<2, double, std::array>({833, 213}, {875, 300});
    }
    if(Layer_name=="TBPS_PosTilted_L1D11"){
      Range = Acts::RangeXD<2, double, std::array>({987, 213}, {1029, 300});
    }
    if(Layer_name=="TBPS_PosTilted_L1D12"){
      Range = Acts::RangeXD<2, double, std::array>({1162, 213}, {1204, 300});
    }
    // TBPS tilted L1 neg
    if(Layer_name=="TBPS_NegTilted_L1D1"){
      Range = Acts::RangeXD<2, double, std::array>({-196, 213}, {-147, 300});
    }
    if(Layer_name=="TBPS_NegTilted_L1D2"){
      Range = Acts::RangeXD<2, double, std::array>({-245, 213}, {-197, 300});
    }
    if(Layer_name=="TBPS_NegTilted_L1D3"){
      Range = Acts::RangeXD<2, double, std::array>({-294, 213}, {-246, 300});
    }
    if(Layer_name=="TBPS_NegTilted_L1D4"){
      Range = Acts::RangeXD<2, double, std::array>({-343, 213}, {-295, 300});
    }
    if(Layer_name=="TBPS_NegTilted_L1D5"){
      Range = Acts::RangeXD<2, double, std::array>({-406, 213}, {-350, 300});
    }
    if(Layer_name=="TBPS_NegTilted_L1D6"){
      Range = Acts::RangeXD<2, double, std::array>({-476, 213}, {-420, 300});
    }
    if(Layer_name=="TBPS_NegTilted_L1D7"){
      Range = Acts::RangeXD<2, double, std::array>({-553, 213}, {-504, 300});
    }
    if(Layer_name=="TBPS_NegTilted_L1D8"){
      Range = Acts::RangeXD<2, double, std::array>({-630, 213}, {-595, 300});
    }
    if(Layer_name=="TBPS_NegTilted_L1D9"){
      Range = Acts::RangeXD<2, double, std::array>({-742, 213}, {-700, 300});
    }
    if(Layer_name=="TBPS_NegTilted_L1D10"){
      Range = Acts::RangeXD<2, double, std::array>({-875, 213}, {-833, 300});
    }
    if(Layer_name=="TBPS_NegTilted_L1D11"){
      Range = Acts::RangeXD<2, double, std::array>({-1029, 213}, {-987, 300});
    }
    if(Layer_name=="TBPS_NegTilted_L1D12"){
      Range = Acts::RangeXD<2, double, std::array>({-1204, 213}, {-1162, 300});
    }
    // TBPS tilted L2 pos
    if(Layer_name=="TBPS_PosTilted_L2D1"){
      Range = Acts::RangeXD<2, double, std::array>({238, 335}, {294, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2D2"){
      Range = Acts::RangeXD<2, double, std::array>({295, 335}, {350, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2D3"){
      Range = Acts::RangeXD<2, double, std::array>({351, 335}, {406, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2D4"){
      Range = Acts::RangeXD<2, double, std::array>({407, 335}, {462, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2D5"){
      Range = Acts::RangeXD<2, double, std::array>({469, 335}, {525, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2D6"){
      Range = Acts::RangeXD<2, double, std::array>({539, 335}, {595, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2D7"){
      Range = Acts::RangeXD<2, double, std::array>({623, 335}, {672, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2D8"){
      Range = Acts::RangeXD<2, double, std::array>({707, 335}, {749, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2D9"){
      Range = Acts::RangeXD<2, double, std::array>({798, 335}, {847, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2D10"){
      Range = Acts::RangeXD<2, double, std::array>({910, 335}, {952, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2D11"){
      Range = Acts::RangeXD<2, double, std::array>({1029, 335}, {1071, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2D12"){
      Range = Acts::RangeXD<2, double, std::array>({1162, 335}, {1204, 415});
    }
    // TBPS tilted L2 neg
    if(Layer_name=="TBPS_NegTilted_L2D1"){
      Range = Acts::RangeXD<2, double, std::array>({-294, 335}, {-238, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2D2"){
      Range = Acts::RangeXD<2, double, std::array>({-350, 335}, {-295, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2D3"){
      Range = Acts::RangeXD<2, double, std::array>({-406, 335}, {-351, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2D4"){
      Range = Acts::RangeXD<2, double, std::array>({-462, 335}, {-407, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2D5"){
      Range = Acts::RangeXD<2, double, std::array>({-525, 335}, {-469, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2D6"){
      Range = Acts::RangeXD<2, double, std::array>({-595, 335}, {-539, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2D7"){
      Range = Acts::RangeXD<2, double, std::array>({-672, 335}, {-623, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2D8"){
      Range = Acts::RangeXD<2, double, std::array>({-749, 335}, {-707, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2D9"){
      Range = Acts::RangeXD<2, double, std::array>({-847, 335}, {-798, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2D10"){
      Range = Acts::RangeXD<2, double, std::array>({-952, 335}, {-910, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2D11"){
      Range = Acts::RangeXD<2, double, std::array>({-1071, 335}, {-1029, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2D12"){
      Range = Acts::RangeXD<2, double, std::array>({-1204, 335}, {-1162, 415});
    }
    // TBPS tilted L3 pos
    if(Layer_name=="TBPS_PosTilted_L3D1"){
      Range = Acts::RangeXD<2, double, std::array>({336, 490}, {385, 566});
    }
    if(Layer_name=="TBPS_PosTilted_L3D2"){
      Range = Acts::RangeXD<2, double, std::array>({386, 490}, {441, 566});
    }
    if(Layer_name=="TBPS_PosTilted_L3D3"){
      Range = Acts::RangeXD<2, double, std::array>({442, 490}, {504, 566});
    }
    if(Layer_name=="TBPS_PosTilted_L3D4"){
      Range = Acts::RangeXD<2, double, std::array>({505, 490}, {567, 566});
    }
    if(Layer_name=="TBPS_PosTilted_L3D5"){
      Range = Acts::RangeXD<2, double, std::array>({568, 490}, {623, 566});
    }
    if(Layer_name=="TBPS_PosTilted_L3D6"){
      Range = Acts::RangeXD<2, double, std::array>({637, 490}, {700, 566});
    }
    if(Layer_name=="TBPS_PosTilted_L3D7"){
      Range = Acts::RangeXD<2, double, std::array>({707, 490}, {756, 566});
    }
    if(Layer_name=="TBPS_PosTilted_L3D8"){
      Range = Acts::RangeXD<2, double, std::array>({784, 490}, {833, 566});
    }
    if(Layer_name=="TBPS_PosTilted_L3D9"){
      Range = Acts::RangeXD<2, double, std::array>({868, 490}, {917, 566});
    }
    if(Layer_name=="TBPS_PosTilted_L3D10"){
      Range = Acts::RangeXD<2, double, std::array>({952, 490}, {1001, 566});
    }
    if(Layer_name=="TBPS_PosTilted_L3D11"){
      Range = Acts::RangeXD<2, double, std::array>({1050, 490}, {1099, 566});
    }
    if(Layer_name=="TBPS_PosTilted_L3D12"){
      Range = Acts::RangeXD<2, double, std::array>({1155, 490}, {1204, 566});
    }
    // TBPS tilted L3 neg
    if(Layer_name=="TBPS_NegTilted_L3D1"){
      Range = Acts::RangeXD<2, double, std::array>({-385, 490}, {-336, 566});
    }
    if(Layer_name=="TBPS_NegTilted_L3D2"){
      Range = Acts::RangeXD<2, double, std::array>({-441, 490}, {-386, 566});
    }
    if(Layer_name=="TBPS_NegTilted_L3D3"){
      Range = Acts::RangeXD<2, double, std::array>({-504, 490}, {-442, 566});
    }
    if(Layer_name=="TBPS_NegTilted_L3D4"){
      Range = Acts::RangeXD<2, double, std::array>({-567, 490}, {-505, 566});
    }
    if(Layer_name=="TBPS_NegTilted_L3D5"){
      Range = Acts::RangeXD<2, double, std::array>({-623, 490}, {-568, 566});
    }
    if(Layer_name=="TBPS_NegTilted_L3D6"){
      Range = Acts::RangeXD<2, double, std::array>({-700, 490}, {-637, 566});
    }
    if(Layer_name=="TBPS_NegTilted_L3D7"){
      Range = Acts::RangeXD<2, double, std::array>({-756, 490}, {-707, 566});
    }
    if(Layer_name=="TBPS_NegTilted_L3D8"){
      Range = Acts::RangeXD<2, double, std::array>({-833, 490}, {-784, 566});
    }
    if(Layer_name=="TBPS_NegTilted_L3D9"){
      Range = Acts::RangeXD<2, double, std::array>({-917, 490}, {-868, 566});
    }
    if(Layer_name=="TBPS_NegTilted_L3D10"){
      Range = Acts::RangeXD<2, double, std::array>({-1001, 490}, {-952, 566});
    }
    if(Layer_name=="TBPS_NegTilted_L3D11"){
      Range = Acts::RangeXD<2, double, std::array>({-1099, 490}, {-1050, 566});
    }
    if(Layer_name=="TBPS_NegTilted_L3D12"){
      Range = Acts::RangeXD<2, double, std::array>({-1204, 490}, {-1155, 566});
    }


    // ======================= Extra Volumes: =======================
    // All TBPS tilted L1 positive
    if(Layer_name=="TBPS_PosTilted_L1_All"){
      Range = Acts::RangeXD<2, double, std::array>({147, 213}, {1204, 300});
    }
    // All TBPS tilted L1 negative
    if(Layer_name=="TBPS_NegTilted_L1_All"){
      Range = Acts::RangeXD<2, double, std::array>({-1204, 213}, {-147, 300});
    }
    // ------
    // All TBPS tilted L2 positive
    if(Layer_name=="TBPS_PosTilted_L2_All"){
      Range = Acts::RangeXD<2, double, std::array>({238, 335}, {1204, 415});
    }
    // All TBPS tilted L2 negative
    if(Layer_name=="TBPS_NegTilted_L2_All"){
      Range = Acts::RangeXD<2, double, std::array>({-1204, 335}, {-238, 415});
    }
    // ------
    // All TBPS tilted L3 positive
    if(Layer_name=="TBPS_PosTilted_L3_All"){
      Range = Acts::RangeXD<2, double, std::array>({336, 490}, {1204, 566});
    }
    // All TBPS tilted L3 negative
    if(Layer_name=="TBPS_NegTilted_L3_All"){
      Range = Acts::RangeXD<2, double, std::array>({-1204, 490}, {-336, 566});
    }

    // ------------------------------------------------------------------
    // TBPS tilted L1 positive two groups:
    if(Layer_name=="TBPS_PosTilted_L1_1to3"){
      Range = Acts::RangeXD<2, double, std::array>({147, 213}, {294, 300});
    }
    if(Layer_name=="TBPS_PosTilted_L1_3to6"){
      Range = Acts::RangeXD<2, double, std::array>({295, 213}, {476, 300});
    }
    // TBPS tilted L1 negative two groups:
    if(Layer_name=="TBPS_NegTilted_L1_1to3"){
      Range = Acts::RangeXD<2, double, std::array>({-294, 213}, {-147, 300});
    }
    if(Layer_name=="TBPS_NegTilted_L1_3to6"){
      Range = Acts::RangeXD<2, double, std::array>({-476, 213}, {-295, 300});
    }

    // TBPS tilted L2 positive three groups:
    if(Layer_name=="TBPS_PosTilted_L2_1to4"){
      Range = Acts::RangeXD<2, double, std::array>({238, 335}, {462, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2_5to8"){
      Range = Acts::RangeXD<2, double, std::array>({469, 335}, {749, 415});
    }
    if(Layer_name=="TBPS_PosTilted_L2_9to12"){
      Range = Acts::RangeXD<2, double, std::array>({798, 335}, {1204, 415});
    }
    // TBPS tilted L2 negative three groups:
    if(Layer_name=="TBPS_NegTilted_L2_1to4"){
      Range = Acts::RangeXD<2, double, std::array>({-462, 335}, {-238, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2_5to8"){
      Range = Acts::RangeXD<2, double, std::array>({-749, 335}, {-469, 415});
    }
    if(Layer_name=="TBPS_NegTilted_L2_9to12"){
      Range = Acts::RangeXD<2, double, std::array>({-1204, 335}, {-798, 415});
    }
 

    // TBPS tilted L3 positive two groups:
    if(Layer_name=="TBPS_PosTilted_L3_1to6"){
      Range = Acts::RangeXD<2, double, std::array>({336, 490}, {700, 566});
    }
    if(Layer_name=="TBPS_PosTilted_L3_6to12"){
      Range = Acts::RangeXD<2, double, std::array>({707, 490}, {1204, 566});
    }
    // TBPS tilted L3 negative two groups:
    if(Layer_name=="TBPS_NegTilted_L3_1to6"){
      Range = Acts::RangeXD<2, double, std::array>({-700, 490}, {-336, 566});
    }
    if(Layer_name=="TBPS_NegTilted_L3_6to12"){
      Range = Acts::RangeXD<2, double, std::array>({-1204, 490}, {-707, 566});
    }


    // 
    
    std::vector<std::shared_ptr<Acts::Surface>> ActiveSurfaces = surfaces.surfaces(Range);

    return ActiveSurfaces;
  }