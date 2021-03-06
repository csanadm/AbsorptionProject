const int NMAT = 201; //number of materials
const int NAI = 6; //number of a_i values for a given material

double aivalues[NMAT][NAI] = 
{{0.02,0.03,0.03,0.03,0.04,0.07},
{0.01,0.01,0.02,0.02,0.02,0.05},
{0.01,0.01,0.01,0.02,0.02,0.02},
{0.05,0.05,0.05,0.08,0.14,0.2},
{0.1,0.2,0.4,0.6,0.5,0.6},
{0.02,0.03,0.03,0.04,0.05,0.07},
{0.01,0.01,0.02,0.02,0.02,0.02},
{0.05,0.04,0.02,0.04,0.05,0.05},
{0.08,0.09,0.12,0.16,0.22,0.24},
{0.02,0.02,0.03,0.04,0.05,0.05},
{0.01,0.01,0.01,0.02,0.02,0.02},
{0.02,0.02,0.02,0.02,0.02,0.02},
{0.02,0.03,0.04,0.05,0.07,0.08},
{0.01,0.01,0.01,0.02,0.02,0.02},
{0.2,0.45,0.6,0.4,0.45,0.4},
{0.04,0.05,0.06,0.08,0.04,0.06},
{0.03,0.03,0.02,0.03,0.04,0.05},
{0.3,0.2,0.15,0.05,0.05,0.05},
{0.3,0.12,0.08,0.06,0.06,0.05},
{0.08,0.11,0.05,0.03,0.02,0.03},
{0.29,0.1,0.05,0.04,0.07,0.09},
{0.31,0.33,0.14,0.1,0.1,0.12},
{0.15,0.1,0.06,0.04,0.04,0.05},
{0.15,0.102469507659596,0.07,0.0529150262212918,0.04,0.05},
{0.08,0.11,0.05,0.03,0.02,0.03},
{0.3,0.12,0.08,0.06,0.06,0.05},
{0.15,0.01,0.06,0.04,0.04,0.05},
{0.3,0.2,0.1,0.07,0.05,0.02},
{0.1,0.06,0.04,0.03,0.02,0.02},
{0.15,0.05,0.03,0.03,0.02,0.02},
{0.5,0.3,0.1,0.05,0.05,0.05},
{0.4,0.35,0.2,0.15,0.05,0.05},
{0.25,0.05,0.04,0.03,0.03,0.02},
{0.28,0.08,0.07,0.07,0.09,0.09},
{0.14,0.1,0.1,0.08,0.1,0.08},
{0.35,0.2,0.15,0.1,0.05,0.05},
{0.4,0.2,0.15,0.1,0.1,0.05},
{0.3,0.2,0.15,0.1,0.1,0.05},
{0.4,0.25,0.15,0.1,0.1,0.05},
{0.31,0.33,0.14,0.1,0.1,0.12},
{0.05,0.25,0.6,0.15,0.05,0.1},
{0.25,0.15,0.1,0.09,0.08,0.07},
{0.25,0.15,0.1,0.09,0.08,0.07},
{0.2,0.62,0.98,0.62,0.21,0.15},
{0.15,0.2,0.1,0.1,0.1,0.1},
{0.03,0.05,0.04,0.03,0.03,0.02},
{0.06,0.1,0.08,0.09,0.07,0.04},
{0.12,0.1,0.08,0.07,0.1,0.08},
{0.3,0.2,0.2,0.1,0.05,0.05},
{0.12,0.04,0.06,0.05,0.05,0.05},
{0.18,0.34,0.42,0.59,0.83,0.68},
{0.19,0.23,0.25,0.3,0.37,0.42},
{0.2,0.25,0.2,0.2,0.15,0.2},
{0.12,0.04,0.06,0.05,0.05,0.05},
{0.09,0.22,0.54,0.76,0.88,0.93},
{0.18,0.56,0.96,1,1,1},
{0.12,0.28,0.55,0.71,0.74,0.83},
{0.17,0.45,0.8,0.89,0.97,0.94},
{0.3,0.69,0.94,1,1,1},
{0.43,0.86,1,1,1,1},
{0.11,0.32,0.56,0.77,0.89,0.91},
{0.27,0.54,0.94,1,0.96,0.96},
{0.28,0.79,1,1,1,1},
{0.46,1,1,1,1,1},
{0.2,0.55,1,1,1,1},
{0.37,0.85,1,1,1,1},
{0.53,0.92,1,1,1,1},
{0.3,0.8,1,1,1,1},
{0.43,0.97,1,1,1,1},
{0.65,1,1,1,1,1},
{0.15,0.6,0.9,0.9,0.9,0.85},
{0.35,0.95,0.98,0.92,0.9,0.85},
{0.11,0.6,0.96,0.94,0.92,0.82},
{0.34,0.95,0.98,0.82,0.87,0.86},
{0.1,0.4,0.8,0.9,0.9,0.9},
{0.4,0.75,0.9,0.8,0.9,0.85},
{0.2,0.374165738677394,0.7,0.793725393319377,0.9,0.8},
{0.3,0.489897948556636,0.8,0.871779788708135,0.95,0.9},
{0.4,0.565685424949238,0.8,0.848528137423857,0.9,0.8},
{0.15,0.4,0.75,0.85,0.8,0.85},
{0.4,0.8,0.95,0.95,0.8,0.85},
{0.1,0.35,0.5,0.55,0.7,0.7},
{0.3,0.55,0.8,0.85,0.75,0.8},
{0.1,0.35,0.55,0.65,0.75,0.8},
{0.2,0.5,0.7,0.8,0.75,0.8},
{0.1,0.25,0.55,0.7,0.8,0.85},
{0.25,0.5,0.85,0.95,0.9,0.9},
{0.2,0.4,0.65,0.55,0.7,0.7},
{0.05,0.15,0.4,0.35,0.2,0.2},
{0.1,0.25,0.55,0.2,0.1,0.15},
{0.05,0.1,0.2,0.55,0.6,0.55},
{0.03,0.05,0.17,0.52,0.5,0.52},
{0.15,0.4,0.65,0.35,0.35,0.3},
{0.01,0.01,0.01,0.01,0.01,0.01},
{0.15,0.45,0.7,0.85,0.95,0.85},
{0.05,0.0707106781186548,0.1,0.1,0.1,0.1},
{0.05,0.1,0.15,0.25,0.3,0.3},
{0.1,0.264575131106459,0.7,0.748331477354788,0.8,0.8},
{0.05,0.0866025403784439,0.15,0.212132034355964,0.3,0.3},
{0.05,0.0707106781186548,0.1,0.122474487139159,0.15,0.15},
{0.3,0.3,0.3,0.3,0.3,0.3},
{0.3,0.212132034355964,0.15,0.122474487139159,0.1,0.1},
{0.2,0.15,0.1,0.05,0.05,0.05},
{0.2,0.15,0.1,0.08,0.04,0.02},
{0.03,0.15,0.5,0.8,0.85,0.8},
{0.02,0.08,0.3,0.6,0.8,0.9},
{0.3,0.3,0.6,0.8,0.75,0.75},
{0.35,0.35,0.4,0.55,0.7,0.7},
{0.1,0.2,0.4,0.489897948556636,0.6,0.6},
{0.1,0.2,0.45,0.8,0.6,0.75},
{0.2,0.4,0.8,0.8,0.8,0.8},
{0.15,0.3,0.6,0.6,0.6,0.7},
{0.3,0.2,0.15,0.05,0.05,0.05},
{0.05,0.15,0.35,0.4,0.5,0.5},
{0.3,0.45,0.65,0.56,0.59,0.71},
{0.05,0.06,0.39,0.63,0.7,0.73},
{0.03,0.03,0.15,0.4,0.5,0.5},
{0.05,0.25,0.4,0.5,0.6,0.5},
{0.05,0.111803398874989,0.25,0.273861278752583,0.3,0.4},
{0.05,0.15,0.35,0.4,0.5,0.5},
{0.1,0.2,0.4,0.447213595499958,0.5,0.6},
{0.11,0.4,0.7,0.74,0.88,0.89},
{0.01,0.01,0.01,0.01,0.02,0.02},
{0.08,0.07,0.06,0.07,0.08,0.08},
{0.27,0.26,0.52,0.43,0.51,0.58},
{0.15,0.11,0.1,0.07,0.06,0.07},
{0.04,0.04,0.07,0.06,0.06,0.07},
{0.2,0.15,0.1,0.1,0.05,0.1},
{0.02,0.02,0.03,0.04,0.04,0.05},
{0.02,0.02,0.04,0.05,0.05,0.1},
{0.01,0.02,0.05,0.15,0.3,0.4},
{0.03,0.09,0.25,0.31,0.33,0.44},
{0.03,0.09,0.2,0.54,0.7,0.72},
{0.08,0.08,0.3,0.6,0.75,0.8},
{0.05,0.05,0.05,0.05,0.05,0.05},
{0.05,0.05,0.1,0.2,0.45,0.65},
{0.5,0.1,0.3,0.5,0.65,0.7},
{0.15,0.25,0.5,0.6,0.7,0.7},
{0.05,0.05,0.1,0.1,0.05,0.05},
{0.1,0.15,0.25,0.3,0.3,0.3},
{0.2,0.25,0.3,0.3,0.3,0.3},
{0.03,0.05,0.05,0.25,0.35,0.5},
{0.02,0.02,0.02,0.0316227766016838,0.05,0.05},
{0.05,0.05,0.15,0.25,0.25,0.25},
{0.05,0.05,0.05,0.1,0.05,0.05},
{0.02,0.04,0.05,0.05,0.1,0.05},
{0.03,0.03,0.03,0.0387298334620742,0.05,0.05},
{0.13,0.09,0.08,0.09,0.11,0.11},
{0.3,0.25,0.15,0.1,0.1,0.07},
{0.14,0.1,0.06,0.08,0.1,0.1},
{0.35,0.39,0.44,0.49,0.54,0.57},
{0.42,0.72,0.83,0.88,0.89,0.8},
{0.06,0.4,0.75,0.95,0.96,0.83},
{0.45,0.7,0.8,0.8,0.65,0.45},
{0.12,0.45,0.87,0.98,1,1},
{0.3,0.7,0.85,0.9,0.7,0.65},
{0.45,0.7,0.88,0.52,0.42,0.35},
{0.5,0.75,0.85,0.65,0.7,0.7},
{0.3,0.4,0.5,0.85,0.5,0.65},
{0.25,0.8,0.85,0.65,0.7,0.75},
{0.15,0.4,0.95,0.6,0.7,0.6},
{0.3,0.2,0.15,0.05,0.05,0.05},
{0.25,0.7,0.85,0.55,0.4,0.3},
{0.4,0.35,0.2,0.15,0.05,0.05},
{0.4,0.2,0.15,0.1,0.1,0.05},
{0.2,0.35,0.55,0.3,0.25,0.3},
{0.4,0.9,0.8,0.5,0.4,0.3},
{0.27,0.87,1,1,0.98,0.96},
{0.5,0.35,0.15,0.05,0.05,0.05},
{0.1,0.3,0.6,0.75,0.8,0.8},
{0.2,0.35,0.65,0.85,0.9,0.8},
{0.1,0.36,0.74,0.91,0.61,0.5},
{0.2,0.22,0.18,0.15,0.15,0.16},
{0.12,0.22,0.37,0.4,0.42,0.37},
{0.28,0.303973683071413,0.33,0.349428104193123,0.37,0.37},
{0.3,0.41,0.49,0.84,0.87,0.84},
{0.33,0.4,0.44,0.45,0.45,0.45},
{0.15,0.38,0.42,0.43,0.45,0.45},
{0.07,0.0989949493661167,0.14,0.14,0.14,0.14},
{0.4,0.5,0.58,0.61,0.58,0.5},
{0.44,0.6,0.77,0.89,0.82,0.7},
{0.49,0.66,0.8,0.88,0.82,0.7},
{0.3,0.346410161513775,0.4,0.414728827066554,0.43,0.4},
{0.16,0.25298221281347,0.4,0.419523539268061,0.44,0.4},
{0.16,0.24,0.56,0.69,0.81,0.78},
{0.24,0.4,0.78,0.98,0.96,0.87},
{0.08,0.109544511501033,0.15,0.16431676725155,0.18,0.2},
{0.07,0.12,0.26,0.42,0.5,0.55},
{0.32,0.62,0.74,0.76,0.81,0.9},
{0.12,0.183303027798234,0.28,0.299332590941915,0.32,0.37},
{0.33,0.51,0.64,0.71,0.77,0.81},
{0.6,0.74,0.88,0.96,0.93,0.85},
{0.13,0.33,0.59,0.58,0.61,0.62},
{0.37,0.48,0.68,0.73,0.77,0.74},
{0.27,0.53,0.67,0.93,0.87,0.8},
{0.37,0.8,1,1,1,1},
{0.2,0.244948974278318,0.3,0.346410161513775,0.4,0.5},
{0.4,0.3,0.2,0.17,0.15,0.1},
{0.5,0.4,0.45,0.45,0.6,0.7},
{0.01,0.01,0.01,0.01,0.02,0.02},
{0.6,0.6,0.6,0.6,0.6,0.6}};

const char *materialnames[NMAT] = 
{"Rough concrete",
"Smooth unpainted concrete",
"Smooth concrete, painted or glazed",
"Porous concrete blocks (no surface finish)",
"Clinker concrete (no surface finish)",
"Smooth brickwork with flush pointing",
"Smooth brickwork with flush pointing, painted",
"Standard brickwork",
"Brickwork, 10mm flush pointing",
"Lime cement plaster on masonry wall",
"Glaze plaster on masonry wall",
"Painted plaster surface on masonry wall",
"Plaster on masonry wall with wall paper on backing paper",
"Ceramic tiles with smooth surface",
"Breeze block",
"Plaster on solid wall",
"Plaster, lime, or gypsum on solid backing",
"Plasterboard on battens, 18mm airspace with glass wool",
"Plasterboard on frame, 100mm airspace",
"Plasterboard on frame, 100mm airspace with glass wool",
"Plasterboard on 50mm battens",
"Plasterboard on 25mm battens",
"2 x plasterboard on frame, 50mm airspace with mineral wool, 2x13 mm",
"Plasterboard on cellular core partition",
"Plasterboard on frame 100mm cavity, 13 mm",
"Plasterboard on frame, 100mm cavity with mineral wool, 13 mm",
" 2 x 13mm plasterboard on steel frame, 50mm mineral wool in cavity, surface painted, 26 mm",
"4mm glass",
"6mm glass",
"Double glazing, 2-3mm glass, 10mm air gap",
"3-4mm plywood, 75mm cavity containing mineral wool",
"5mm plywood on battens, 50mm airspace filled",
"12mm plywood over 50mm airgap",
"12mm plywood over 150mm airgap",
"12mm plywood over 200mm airgap containing 50mm mineral wool",
"12mm plywood in framework with 30mm airspace behind",
"12mm plywood in framework with 30mm airspace containing glass wool",
"Plywood, hardwood panels over 25mm airspace on solid backing",
"Plywood, hardwood panels over 25mm airspace on solid backing with absorbent material in air space",
"12mm wood panelling on 25mm battens",
"Timber boards, 100mm wide, 10mm gaps, 500mm airspace with mineral wool, 22 mm",
"t & g board on frame, 50mm airspace with mineral wool, 16 mm",
"16-22mm t&g wood on 50mm cavity filled with mineral wool",
"Cedar, slotted and profiled on battens mineral wool in airspace",
"Wood boards on on joists or battens",
"20mm dense veneered chipboard over 100mm airgap",
"20mm dense veneered chipboard over 200mm airgap",
"20mm dense veneered chipboard over 250mm airgap containing 50mm mineral wool",
"6mm wood fibre board, cavity > 100mm, empty",
"22mm chipboard, 50mm cavity filled with mineral wool",
"Acoustic timber wall panelling",
"Hardwood, mahogany",
"Chipboard on 16mm battens, 20 mm thick",
"Chipboard on frame, 50mm airspace with mineral wool, 22 mm thick",
"Melamine based foam 25mm",
"Melamine based foam 50mm",
"Glass wool 25mm 16 kg/m3",
"Glass wool 50mm, 16 kg/m3",
"Glass wool 75mm, 16 kg/m3",
"Glass wool 100mm, 16 kg/m3",
"Glass wool 25mm, 24 kg/m3",
"Glass wool 50mm, 24 kg/m3",
"Glass wool 75mm, 24 kg/m3",
"Glass wool 100mm, 24 kg/m3",
"Glass wool 50mm, 33 kg/m3",
"Glass wool 75mm, 33 kg/m3",
"Glass wool 100mm, 33 kg/m3",
"Glass wool 50mm, 48 kg/m3",
"Glass wool 75mm, 48 kg/m3",
"Glass wool 100mm, 48 kg/m3",
"Rock wool 50mm, 33 kg/m3 direct to masonry",
"Rock wool 100mm, 33 kg/m3 direct to masonry",
"Rock wool 50mm, 60 kg/m3 direct to masonry",
"Rock wool 75mm, 60 kg/m3 direct to masonry",
"Rock wool 30mm, 100 kg/m3 direct to masonry",
"Rock wool 30mm, 200 kg/m3 over 300mm air gap",
"Glass wool or mineral wool on solid backing, 25 mm",
"Glass wool or mineral wool on solid backing, 50 mm",
"Glass wool or mineral wool over air space on solid backing, 25 mm",
"Fibreglass super fine mat, 50 mm",
"Fibreglass scrim-covered sewn sheet, 40 mm",
"Fibreglass bitumen bonded mat, 25 mm",
"Fibreglass bitumen bonded mat, 50 mm",
"Fibreglass resin-bonded mat, 25 mm",
"Fibreglass resin-bonded mat, 50 mm",
"Fibreglass resin-bonded board, 25 mm",
"Flexible polyurethane foam 50mm",
"Rigid polyurethane foam 50mm",
"12mm expanded polystyrene on 45mm battens",
"25mm expanded polystyrene on 50mm battens",
"Cork tiles 25mm on solid backing",
"Cork board, 25mm on solid backing",
"Cork board, 25mm, 2.9kg/m2, on battens",
"Glass blocks or glazed tiles as wall finish",
"Muslin covered cotton felt , 25 mm",
"Pin up boarding- medium hardboard on solid backing",
"Fibreboard on solid backing , 12 mm",
"25mm thick hair felt, covered by scrim cloth on solid backing",
"Fibreboard on solid backing, soft 12 mm",
"Fibreboard on solid backing - painted",
"Fibreboard over airspace on solid wall, 12 mm",
"Fibreboard over airspace on solid wall - painted",
"Plaster on lath, deep air space",
"Plaster decorative panels, walls",
"Acoustic plaster to solid backing, 25 mm",
"9mm acoustic plaster to solid backing",
"9mm acoustic plaster on plasterboard, 75mm airspace ",
"12.5mm acoustic plaster on plaster backing over 75mm air space",
"Woodwool slabs, unplastered on solid backing, 25 mm",
"Woodwool slabs, unplastered on solid backing, 50 mm",
"Woodwool slabs, unplastered on solid backing, 75 mm",
"Woodwool slabs, unplastered over 20mm airspace on solid backing, 25 mm thick",
"Plasterboard backed with 25mm thick bitumen-bonded fibreglass on 50mm battens, 10 mm thick",
"Curtains hung in folds against soild wall",
"Cotton Curtains (0.5kg/m2),draped to 75% area approx. 130mm from wall",
"Lightweight curtains (0.2 kg/m2) hung 90mm from wall",
"Curtains of close-woven glass mat hung 50mm from wall",
"Curtains, medium velour, 50% gather, over soild backing",
"Curtains (medium fabrics) hung straight and close to wall",
"Curtains in folds against wall",
"Curtains (medium fabrics) double widths in folds spaced away from wall",
"Acoustic banner, 0.5 kg/m2 wool serge, 100mm from wall",
"Smooth marble or terrazzo slabs",
"Raised computer floor, steel-faced 45mm chipboard 800mm above concrete floor, no carpet",
"Raised computer floor, steel-faced 45mm chipboard 800mm above concrete floor, office-grade carpet tiles",
"Wooden floor on joists",
"Parquet fixed in asphalt, on concrete",
"Parquet on counterfloor",
"Linoleum or vinyl stuck to concrete",
"Layer of rubber, cork, linoleum + underlay, or vinyl+underlay stuck to concrete",
"5mm needle-felt stuck to concrete",
"6mm pile carpet bonded to closed-cell foam underlay",
"6mm pile carpet bonded to open-cell foam underlay",
"9mm pile carpet, tufted on felt underlay",
"Composition flooring",
"Haircord carpet on felt underlay, 6 mm",
"Medium pile carpet on sponge rubber underlay, 10 mm",
"Thick pile carpet on sponge rubber underlay, 15 mm",
"Rubber floor tiles, 6 mm",
"Carpet, thin, over thin felt on concrete",
"Carpet, thin, over thin felt on wood floor",
"Carpet, needlepunch, 5 mm",
"Stone floor, plain or tooled or granolithic finish",
"Corkfloor tiles, 14 mm",
"Sheet rubber (hard) , 6 mm",
"Woodblock/linoleum/rubber/cork tiles (thin) on solid floor (or wall)",
"Floor tiles, plastic, or linoleum",
"Steel decking",
"Wood hollowcore door",
"Solid timber door",
"Acoustic door, steel frame, double seals, absorbant in airspace, Double sheet steel skin.",
"Mineral wool tiles, 180mm airspace",
"Mineral wool tiles, glued/screwed to soffit",
"Gypsum plaster tiles, 17% perforated, 22mm",
"Metal ceiling, 32.5% perforated, backed by 30mm rockwool",
"Perforated underside of structural steel decking (typical, depends on perforations)",
"12% perforated plaster tiles, absorbent felt glued to back, 200mm ceiling void",
"100mm woodwool slabs on 25mm cavity, pre-screeded surface facing cavity",
"50mm woodwool slabs on 25mm cavity, pre-screeded surface facing cavity",
"100mm woodwool fixed directly to concrete, pre-screeded surface facing backing",
"75mm woodwool fixed directly to concrete, pre-screeded surface facing backing",
"Plasterboard 10mm thick backed with 25mm thick bitumen",
"Plasterboard 10mm thick, perforated 8mm diameter holes 2755m2 14% open area backed with 25mm thick bitumen- bonded fibreglass on 90mm battens",
"Plywood, 5mm, on battens 50mm airspace filled with glass wool",
"Plywood, 12mm, with 30mm thick fibreglass backing between 30mm battens",
"Plywood 12mm thick perforated 5mm diameter holes 6200 m2 11% open area with 60mm deep air space behind",
"Plywood 12mm thick perforated 5mm diameter holes 6200 m2 11% open area backed with 60mm thick fibreglass between mounting battens",
"Hardboard, 25% perforated over 50mm mineral wool",
"0.8mm unperforated metal panels backed with 25mm thick resin bonded fibreglass, mounted on  22mm diameter pipes 135mm from wall.",
"0.8mm perforated metal tiles 2mm diameter holes 29440/m2. 13% open area backed with 25mm thick resin-bonded fibreglass slab. No airspace.",
"50mm mineral wool (96 kg/m3) behind 25% open area perforated steel.",
"Wood panels, 18mm alternate 15mm slot & 35mm wooden slat",
"25mm rockwool backing, 32mm airspace behind plaster decorative panels, ceilings",
"Children, standing (per child) in m2 units",
"Children, seated in plastic or metal chairs (per child)in m2 units",
"Students seated in tablet arm chairs",
"Adults per person seated",
"Adults per person standing",
"Empty plastic or metal chairs (per chair) in m2 units",
"Seats, leather covers, per m2",
"Cloth-upholstered seats, per m2",
"Floor and cloth-upholstered seats, per m2",
"Adults in plastic and metal chairs in m2 units",
"Adults in wooden or padded chairs or seats (per item) in m2",
"Adults on timber seats, 1 per m2 per item",
"Adults on timber seats, 2 per m2 per item",
"Wooden or padded chairs or seats (per item) in m2",
"Seating, slighty upholstered, unoccupied",
"Seating, slighty upholstered, occupied",
"Fully upholstered seats (per item) in m2",
"Upholstered tip-up theatre seats, empty ",
"Areas with audience, orchestra, or seats, including narrow aisles",
"Auditorium seat, unoccupied",
"Auditorium seat, occupied",
"Orchestra with instruments on podium, 1.5 m2 per person",
"Orchestral player with instrument (average) per person",
"Prosenium opening with average stage set per m2of opening",
"Wood platform with large space beneath",
"Adult office furniture per desk",
"Water surface, ie swimming pool",
"Ventilation grille per m2"};

double frequenciesOtherMat[NAI] = {125,250,500,1000,2000,4000};

double AlphaOtherMaterial(const double *x, const double *par)
{
  int nmat = int(par[0]);
  if(nmat<0 || nmat>=NAI) return -9999;
  double f = x[0];
  
  int range = -1;
  for(int ifr=0; ifr<NAI; ifr++)
    if(f<frequenciesOtherMat[ifr])
    {
      range=ifr;
      break;
    }
  if(range==-1) range=NAI; //this means f larger than largest frequency point
  if(range==0) range=1; //first range is up to 250
  if(range>NAI-2) range=NAI-1; //last range is from 2000
    
  double n = log(aivalues[nmat][range]/aivalues[nmat][range-1]) / log(frequenciesOtherMat[range]/frequenciesOtherMat[range-1]);
  double c = aivalues[nmat][range]/pow(frequenciesOtherMat[range],n);
  
  return c*pow(f,n);
  
}