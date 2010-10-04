#//first offpoint flatfield
#//write offpoint

#//04.08.2010

#recreate_flatfield.csh
#See mail by Richard Wachter 09/29/2010 16:49 "recreate flatfield series"

#!/bin/csh

write_offpoint pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_01.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:07:40.00_UTC" fsn_list_offpoint="3614678,3614712,3614746,3614882,3614916,3614950,3615052,3615086,3615120,3615222,3615256,3615290,3615460,3615494,3615528,3615698,3615732,3615766,3615902,3615936,3615970,3616072,3616106,3616140,3616276,3616310,3616344,3616548,3616582,3616616,3616854,3616888,3616922,3617092,3617126,3617160,3617296,3617330,3617364,3617534,3617568,3617602,3617772,3617806,3617840,3618078,3618112,3618146,3618316,3618350,3618384,3618588,3618622,3618656,3618860,3618894,3618928,3619132,3619166,3619200" nx=4096 ny=4096 focus=01

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_02.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:07:43.00_UTC" fsn_list_offpoint="3614680,3614714,3614748,3614884,3614918,3614952,3615054,3615088,3615122,3615224,3615258,3615292,3615462,3615496,3615530,3615700,3615734,3615768,3615904,3615938,3615972,3616074,3616108,3616142,3616278,3616312,3616346,3616550,3616584,3616618,3616856,3616890,3616924,3617094,3617128,3617162,3617298,3617332,3617366,3617536,3617570,3617604,3617774,3617808,3617842,3618080,3618114,3618148,3618318,3618352,3618386,3618590,3618624,3618658,3618862,3618896,3618930,3619134,3619168,3619202" nx=4096 ny=4096 focus=02

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_03.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:07:47.00_UTC" fsn_list_offpoint="3614682,3614716,3614750,3614886,3614920,3614954,3615056,3615090,3615124,3615226,3615260,3615294,3615464,3615498,3615532,3615702,3615736,3615770,3615906,3615940,3615974,3616076,3616110,3616144,3616280,3616314,3616348,3616552,3616586,3616620,3616858,3616892,3616926,3617096,3617130,3617164,3617300,3617334,3617368,3617538,3617572,3617606,3617776,3617810,3617844,3618082,3618116,3618150,3618320,3618354,3618388,3618592,3618626,3618660,3618864,3618898,3618932,3619136,3619170,3619204" nx=4096 ny=4096 focus=03

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_04.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:07:51.00_UTC" fsn_list_offpoint="3614684,3614718,3614752,3614888,3614922,3614956,3615058,3615092,3615126,3615228,3615262,3615296,3615466,3615500,3615534,3615704,3615738,3615772,3615908,3615942,3615976,3616078,3616112,3616146,3616282,3616316,3616350,3616554,3616588,3616622,3616860,3616894,3616928,3617098,3617132,3617166,3617302,3617336,3617370,3617540,3617574,3617608,3617778,3617812,3617846,3618084,3618118,3618152,3618322,3618356,3618390,3618594,3618628,3618662,3618866,3618900,3618934,3619138,3619172,3619206" nx=4096 ny=4096 focus=04

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_05.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:07:55.00_UTC" fsn_list_offpoint="3614686,3614720,3614754,3614890,3614924,3614958,3615060,3615094,3615128,3615230,3615264,3615298,3615468,3615502,3615536,3615706,3615740,3615774,3615910,3615944,3615978,3616080,3616114,3616148,3616284,3616318,3616352,3616556,3616590,3616624,3616862,3616896,3616930,3617100,3617134,3617168,3617304,3617338,3617372,3617542,3617576,3617610,3617780,3617814,3617848,3618086,3618120,3618154,3618324,3618358,3618392,3618596,3618630,3618664,3618868,3618902,3618936,3619140,3619174,3619208" nx=4096 ny=4096 focus=05

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_06.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:07:58.00_UTC" fsn_list_offpoint="3614688,3614722,3614756,3614892,3614926,3614960,3615062,3615096,3615130,3615232,3615266,3615300,3615470,3615504,3615538,3615708,3615742,3615776,3615912,3615946,3615980,3616082,3616116,3616150,3616286,3616320,3616354,3616558,3616592,3616626,3616864,3616898,3616932,3617102,3617136,3617170,3617306,3617340,3617374,3617544,3617578,3617612,3617782,3617816,3617850,3618088,3618122,3618156,3618326,3618360,3618394,3618598,3618632,3618666,3618870,3618904,3618938,3619142,3619176,3619210" nx=4096 ny=4096 focus=06

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_07.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:08:02.00_UTC" fsn_list_offpoint="3614690,3614724,3614758,3614894,3614928,3614962,3615064,3615098,3615132,3615234,3615268,3615302,3615472,3615506,3615540,3615710,3615744,3615778,3615914,3615948,3615982,3616084,3616118,3616152,3616288,3616322,3616356,3616560,3616594,3616628,3616866,3616900,3616934,3617104,3617138,3617172,3617308,3617342,3617376,3617546,3617580,3617614,3617784,3617818,3617852,3618090,3618124,3618158,3618328,3618362,3618396,3618600,3618634,3618668,3618872,3618906,3618940,3619144,3619178,3619212" nx=4096 ny=4096 focus=07

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_08.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:08:06.00_UTC" fsn_list_offpoint="3614692,3614726,3614760,3614896,3614930,3614964,3615066,3615100,3615134,3615236,3615270,3615304,3615474,3615508,3615542,3615712,3615746,3615780,3615916,3615950,3615984,3616086,3616120,3616154,3616290,3616324,3616358,3616562,3616596,3616630,3616868,3616902,3616936,3617106,3617140,3617174,3617310,3617344,3617378,3617548,3617582,3617616,3617786,3617820,3617854,3618092,3618126,3618160,3618330,3618364,3618398,3618602,3618636,3618670,3618874,3618908,3618942,3619146,3619180,3619214" nx=4096 ny=4096 focus=08

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_09.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:08:10.00_UTC" fsn_list_offpoint="3614694,3614728,3614762,3614898,3614932,3614966,3615068,3615102,3615136,3615238,3615272,3615306,3615476,3615510,3615544,3615714,3615748,3615782,3615918,3615952,3615986,3616088,3616122,3616156,3616292,3616326,3616360,3616564,3616598,3616632,3616870,3616904,3616938,3617108,3617142,3617176,3617312,3617346,3617380,3617550,3617584,3617618,3617788,3617822,3617856,3618094,3618128,3618162,3618332,3618366,3618400,3618604,3618638,3618672,3618876,3618910,3618944,3619148,3619182,3619216" nx=4096 ny=4096 focus=09

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_10.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:08:13.00_UTC" fsn_list_offpoint="3614696,3614730,3614764,3614900,3614934,3614968,3615070,3615104,3615138,3615240,3615274,3615308,3615478,3615512,3615546,3615716,3615750,3615784,3615920,3615954,3615988,3616090,3616124,3616158,3616294,3616328,3616362,3616566,3616600,3616634,3616872,3616906,3616940,3617110,3617144,3617178,3617314,3617348,3617382,3617552,3617586,3617620,3617790,3617824,3617858,3618096,3618130,3618164,3618334,3618368,3618402,3618606,3618640,3618674,3618878,3618912,3618946,3619150,3619184,3619218" nx=4096 ny=4096 focus=10

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_11.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:08:17.00_UTC" fsn_list_offpoint="3614698,3614732,3614766,3614902,3614936,3614970,3615072,3615106,3615140,3615242,3615276,3615310,3615480,3615514,3615548,3615718,3615752,3615786,3615922,3615956,3615990,3616092,3616126,3616160,3616296,3616330,3616364,3616568,3616602,3616636,3616874,3616908,3616942,3617112,3617146,3617180,3617316,3617350,3617384,3617554,3617588,3617622,3617792,3617826,3617860,3618098,3618132,3618166,3618336,3618370,3618404,3618608,3618642,3618676,3618880,3618914,3618948,3619152,3619186,3619220" nx=4096 ny=4096 focus=11

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_12.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:08:21.00_UTC" fsn_list_offpoint="3614700,3614734,3614768,3614904,3614938,3614972,3615074,3615108,3615142,3615244,3615278,3615312,3615482,3615516,3615550,3615720,3615754,3615788,3615924,3615958,3615992,3616094,3616128,3616162,3616298,3616332,3616366,3616570,3616604,3616638,3616876,3616910,3616944,3617114,3617148,3617182,3617318,3617352,3617386,3617556,3617590,3617624,3617794,3617828,3617862,3618100,3618134,3618168,3618338,3618372,3618406,3618610,3618644,3618678,3618882,3618916,3618950,3619154,3619188,3619222" nx=4096 ny=4096 focus=12

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_13.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:08:25.00_UTC" fsn_list_offpoint="3614702,3614736,3614770,3614906,3614940,3614974,3615076,3615110,3615144,3615246,3615280,3615314,3615484,3615518,3615552,3615722,3615756,3615790,3615926,3615960,3615994,3616096,3616130,3616164,3616300,3616334,3616368,3616572,3616606,3616640,3616878,3616912,3616946,3617116,3617150,3617184,3617320,3617354,3617388,3617558,3617592,3617626,3617796,3617830,3617864,3618102,3618136,3618170,3618340,3618374,3618408,3618612,3618646,3618680,3618884,3618918,3618952,3619156,3619190,3619224" nx=4096 ny=4096 focus=13

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_14.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:08:28.00_UTC" fsn_list_offpoint="3614704,3614738,3614772,3614908,3614942,3614976,3615078,3615112,3615146,3615248,3615282,3615316,3615486,3615520,3615554,3615724,3615758,3615792,3615928,3615962,3615996,3616098,3616132,3616166,3616302,3616336,3616370,3616574,3616608,3616642,3616880,3616914,3616948,3617118,3617152,3617186,3617322,3617356,3617390,3617560,3617594,3617628,3617798,3617832,3617866,3618104,3618138,3618172,3618342,3618376,3618410,3618614,3618648,3618682,3618886,3618920,3618954,3619158,3619192,3619226" nx=4096 ny=4096 focus=14

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_15.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:08:32.00_UTC" fsn_list_offpoint="3614706,3614740,3614774,3614910,3614944,3614978,3615080,3615114,3615148,3615250,3615284,3615318,3615488,3615522,3615556,3615726,3615760,3615794,3615930,3615964,3615998,3616100,3616134,3616168,3616304,3616338,3616372,3616576,3616610,3616644,3616882,3616916,3616950,3617120,3617154,3617188,3617324,3617358,3617392,3617562,3617596,3617630,3617800,3617834,3617868,3618106,3618140,3618174,3618344,3618378,3618412,3618616,3618650,3618684,3618888,3618922,3618956,3619160,3619194,3619228" nx=4096 ny=4096 focus=15

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_side_16.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.08_22:08:36.00_UTC" fsn_list_offpoint="3614708,3614742,3614776,3614912,3614946,3614980,3615082,3615116,3615150,3615252,3615286,3615320,3615490,3615524,3615558,3615728,3615762,3615796,3615932,3615966,3616000,3616102,3616136,3616170,3616306,3616340,3616374,3616578,3616612,3616646,3616884,3616918,3616952,3617122,3617156,3617190,3617326,3617360,3617394,3617564,3617598,3617632,3617802,3617836,3617870,3618108,3618142,3618176,3618346,3618380,3618414,3618618,3618652,3618686,3618890,3618924,3618958,3619162,3619196,3619230" nx=4096 ny=4096 focus=16

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_01.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:07:38.00_UTC" fsn_list_offpoint="3614677,3614711,3614745,3614881,3614915,3614949,3615051,3615085,3615119,3615221,3615255,3615289,3615459,3615493,3615527,3615697,3615731,3615765,3615901,3615935,3615969,3616071,3616105,3616139,3616275,3616309,3616343,3616547,3616581,3616615,3616853,3616887,3616921,3617091,3617125,3617159,3617295,3617329,3617363,3617533,3617567,3617601,3617771,3617805,3617839,3618077,3618111,3618145,3618315,3618349,3618383,3618587,3618621,3618655,3618859,3618893,3618927,3619131,3619165,3619199" nx=4096 ny=4096 focus=01

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_02.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:07:42.00_UTC" fsn_list_offpoint="3614679,3614713,3614747,3614883,3614917,3614951,3615053,3615087,3615121,3615223,3615257,3615291,3615461,3615495,3615529,3615699,3615733,3615767,3615903,3615937,3615971,3616073,3616107,3616141,3616277,3616311,3616345,3616549,3616583,3616617,3616855,3616889,3616923,3617093,3617127,3617161,3617297,3617331,3617365,3617535,3617569,3617603,3617773,3617807,3617841,3618079,3618113,3618147,3618317,3618351,3618385,3618589,3618623,3618657,3618861,3618895,3618929,3619133,3619167,3619201" nx=4096 ny=4096 focus=02

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_03.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:07:45.00_UTC" fsn_list_offpoint="3614681,3614715,3614749,3614885,3614919,3614953,3615055,3615089,3615123,3615225,3615259,3615293,3615463,3615497,3615531,3615701,3615735,3615769,3615905,3615939,3615973,3616075,3616109,3616143,3616279,3616313,3616347,3616551,3616585,3616619,3616857,3616891,3616925,3617095,3617129,3617163,3617299,3617333,3617367,3617537,3617571,3617605,3617775,3617809,3617843,3618081,3618115,3618149,3618319,3618353,3618387,3618591,3618625,3618659,3618863,3618897,3618931,3619135,3619169,3619203" nx=4096 ny=4096 focus=03

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_04.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:07:49.00_UTC" fsn_list_offpoint="3614683,3614717,3614751,3614887,3614921,3614955,3615057,3615091,3615125,3615227,3615261,3615295,3615465,3615499,3615533,3615703,3615737,3615771,3615907,3615941,3615975,3616077,3616111,3616145,3616281,3616315,3616349,3616553,3616587,3616621,3616859,3616893,3616927,3617097,3617131,3617165,3617301,3617335,3617369,3617539,3617573,3617607,3617777,3617811,3617845,3618083,3618117,3618151,3618321,3618355,3618389,3618593,3618627,3618661,3618865,3618899,3618933,3619137,3619171,3619205" nx=4096 ny=4096 focus=04

write_offpoint   pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_05.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:07:53.00_UTC" fsn_list_offpoint="3614685,3614719,3614753,3614889,3614923,3614957,3615059,3615093,3615127,3615229,3615263,3615297,3615467,3615501,3615535,3615705,3615739,3615773,3615909,3615943,3615977,3616079,3616113,3616147,3616283,3616317,3616351,3616555,3616589,3616623,3616861,3616895,3616929,3617099,3617133,3617167,3617303,3617337,3617371,3617541,3617575,3617609,3617779,3617813,3617847,3618085,3618119,3618153,3618323,3618357,3618391,3618595,3618629,3618663,3618867,3618901,3618935,3619139,3619173,3619207" nx=4096 ny=4096 focus=05

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_06.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:07:57.00_UTC" fsn_list_offpoint="3614687,3614721,3614755,3614891,3614925,3614959,3615061,3615095,3615129,3615231,3615265,3615299,3615469,3615503,3615537,3615707,3615741,3615775,3615911,3615945,3615979,3616081,3616115,3616149,3616285,3616319,3616353,3616557,3616591,3616625,3616863,3616897,3616931,3617101,3617135,3617169,3617305,3617339,3617373,3617543,3617577,3617611,3617781,3617815,3617849,3618087,3618121,3618155,3618325,3618359,3618393,3618597,3618631,3618665,3618869,3618903,3618937,3619141,3619175,3619209" nx=4096 ny=4096 focus=06

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_07.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:08:00.00_UTC" fsn_list_offpoint="3614689,3614723,3614757,3614893,3614927,3614961,3615063,3615097,3615131,3615233,3615267,3615301,3615471,3615505,3615539,3615709,3615743,3615777,3615913,3615947,3615981,3616083,3616117,3616151,3616287,3616321,3616355,3616559,3616593,3616627,3616865,3616899,3616933,3617103,3617137,3617171,3617307,3617341,3617375,3617545,3617579,3617613,3617783,3617817,3617851,3618089,3618123,3618157,3618327,3618361,3618395,3618599,3618633,3618667,3618871,3618905,3618939,3619143,3619177,3619211" nx=4096 ny=4096 focus=07

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_08.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:08:04.00_UTC" fsn_list_offpoint="3614691,3614725,3614759,3614895,3614929,3614963,3615065,3615099,3615133,3615235,3615269,3615303,3615473,3615507,3615541,3615711,3615745,3615779,3615915,3615949,3615983,3616085,3616119,3616153,3616289,3616323,3616357,3616561,3616595,3616629,3616867,3616901,3616935,3617105,3617139,3617173,3617309,3617343,3617377,3617547,3617581,3617615,3617785,3617819,3617853,3618091,3618125,3618159,3618329,3618363,3618397,3618601,3618635,3618669,3618873,3618907,3618941,3619145,3619179,3619213" nx=4096 ny=4096 focus=08

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_09.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:08:08.00_UTC" fsn_list_offpoint="3614693,3614727,3614761,3614897,3614931,3614965,3615067,3615101,3615135,3615237,3615271,3615305,3615475,3615509,3615543,3615713,3615747,3615781,3615917,3615951,3615985,3616087,3616121,3616155,3616291,3616325,3616359,3616563,3616597,3616631,3616869,3616903,3616937,3617107,3617141,3617175,3617311,3617345,3617379,3617549,3617583,3617617,3617787,3617821,3617855,3618093,3618127,3618161,3618331,3618365,3618399,3618603,3618637,3618671,3618875,3618909,3618943,3619147,3619181,3619215" nx=4096 ny=4096 focus=09

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_10.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:08:12.00_UTC" fsn_list_offpoint="3614695,3614729,3614763,3614899,3614933,3614967,3615069,3615103,3615137,3615239,3615273,3615307,3615477,3615511,3615545,3615715,3615749,3615783,3615919,3615953,3615987,3616089,3616123,3616157,3616293,3616327,3616361,3616565,3616599,3616633,3616871,3616905,3616939,3617109,3617143,3617177,3617313,3617347,3617381,3617551,3617585,3617619,3617789,3617823,3617857,3618095,3618129,3618163,3618333,3618367,3618401,3618605,3618639,3618673,3618877,3618911,3618945,3619149,3619183,3619217" nx=4096 ny=4096 focus=10

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_11.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:08:15.00_UTC" fsn_list_offpoint="3614697,3614731,3614765,3614901,3614935,3614969,3615071,3615105,3615139,3615241,3615275,3615309,3615479,3615513,3615547,3615717,3615751,3615785,3615921,3615955,3615989,3616091,3616125,3616159,3616295,3616329,3616363,3616567,3616601,3616635,3616873,3616907,3616941,3617111,3617145,3617179,3617315,3617349,3617383,3617553,3617587,3617621,3617791,3617825,3617859,3618097,3618131,3618165,3618335,3618369,3618403,3618607,3618641,3618675,3618879,3618913,3618947,3619151,3619185,3619219" nx=4096 ny=4096 focus=11

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_12.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:08:19.00_UTC" fsn_list_offpoint="3614699,3614733,3614767,3614903,3614937,3614971,3615073,3615107,3615141,3615243,3615277,3615311,3615481,3615515,3615549,3615719,3615753,3615787,3615923,3615957,3615991,3616093,3616127,3616161,3616297,3616331,3616365,3616569,3616603,3616637,3616875,3616909,3616943,3617113,3617147,3617181,3617317,3617351,3617385,3617555,3617589,3617623,3617793,3617827,3617861,3618099,3618133,3618167,3618337,3618371,3618405,3618609,3618643,3618677,3618881,3618915,3618949,3619153,3619187,3619221" nx=4096 ny=4096 focus=12

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_13.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:08:23.00_UTC" fsn_list_offpoint="3614701,3614735,3614769,3614905,3614939,3614973,3615075,3615109,3615143,3615245,3615279,3615313,3615483,3615517,3615551,3615721,3615755,3615789,3615925,3615959,3615993,3616095,3616129,3616163,3616299,3616333,3616367,3616571,3616605,3616639,3616877,3616911,3616945,3617115,3617149,3617183,3617319,3617353,3617387,3617557,3617591,3617625,3617795,3617829,3617863,3618101,3618135,3618169,3618339,3618373,3618407,3618611,3618645,3618679,3618883,3618917,3618951,3619155,3619189,3619223" nx=4096 ny=4096 focus=13

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_14.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:08:27.00_UTC" fsn_list_offpoint="3614703,3614737,3614771,3614907,3614941,3614975,3615077,3615111,3615145,3615247,3615281,3615315,3615485,3615519,3615553,3615723,3615757,3615791,3615927,3615961,3615995,3616097,3616131,3616165,3616301,3616335,3616369,3616573,3616607,3616641,3616879,3616913,3616947,3617117,3617151,3617185,3617321,3617355,3617389,3617559,3617593,3617627,3617797,3617831,3617865,3618103,3618137,3618171,3618341,3618375,3618409,3618613,3618647,3618681,3618885,3618919,3618953,3619157,3619191,3619225" nx=4096 ny=4096 focus=14

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_15.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:08:30.00_UTC" fsn_list_offpoint="3614705,3614739,3614773,3614909,3614943,3614977,3615079,3615113,3615147,3615249,3615283,3615317,3615487,3615521,3615555,3615725,3615759,3615793,3615929,3615963,3615997,3616099,3616133,3616167,3616303,3616337,3616371,3616575,3616609,3616643,3616881,3616915,3616949,3617119,3617153,3617187,3617323,3617357,3617391,3617561,3617595,3617629,3617799,3617833,3617867,3618105,3618139,3618173,3618343,3618377,3618411,3618615,3618649,3618683,3618887,3618921,3618955,3619159,3619193,3619227" nx=4096 ny=4096 focus=15

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/tmp20/richard/hmi/flat_front_16.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.08_22:08:34.00_UTC" fsn_list_offpoint="3614707,3614741,3614775,3614911,3614945,3614979,3615081,3615115,3615149,3615251,3615285,3615319,3615489,3615523,3615557,3615727,3615761,3615795,3615931,3615965,3615999,3616101,3616135,3616169,3616305,3616339,3616373,3616577,3616611,3616645,3616883,3616917,3616951,3617121,3617155,3617189,3617325,3617359,3617393,3617563,3617597,3617631,3617801,3617835,3617869,3618107,3618141,3618175,3618345,3618379,3618413,3618617,3618651,3618685,3618889,3618923,3618957,3619161,3619195,3619229" nx=4096 ny=4096 focus=16


#//write dark
write_dark instrument="HMI" file_dark="/tmp20/richard/hmi/dark_front_inorbit1.bin" series_dark="hmi.dark" camera=2 t_obs="2010.04.07_14:37:41.00_UTC" fsn_list_dark=3563653,3563687,3563721,3563755,3563789 nx=4096 ny=4096

write_dark instrument="HMI" file_dark="/tmp20/richard/hmi/dark_side_inorbit1.bin" series_dark="hmi.dark" camera=1 t_obs="2010.04.07_14:37:43.00_UTC" fsn_list_dark=3563654,3563688,3563722,3563756,3563790 nx=4096 ny=4096


#//write_badpix
write_badpix instrument="HMI" file_badpix="/tmp20/richard/hmi/bad_pix_front_inorbit1.bin" series_badpix="hmi.bad_pixel_list"  camera=2 t_obs="2010.04.08_22:07:38.00_UTC" nbad=29

write_badpix instrument="HMI" file_badpix="/tmp20/richard/hmi/bad_pix_side_inorbit1.bin" series_badpix="hmi.bad_pixel_list" camera=1 t_obs="2010.04.08_22:07:40.00_UTC" nbad=45

#//write_flatfield
#//write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:07:40.00_UTC" camera=1 focus=01 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:07:40.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_01.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:07:43.00_UTC" camera=1 focus=02 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:07:43.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_02.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:07:47.00_UTC" camera=1 focus=03 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:07:47.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_03.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:07:51.00_UTC" camera=1 focus=04 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:07:51.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_04.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:07:55.00_UTC" camera=1 focus=05 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:07:55.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_05.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:07:58.00_UTC" camera=1 focus=06 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:07:58.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_06.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:02.00_UTC" camera=1 focus=07 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:08:02.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_07.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:06.00_UTC" camera=1 focus=08 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:08:06.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_08.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:10.00_UTC" camera=1 focus=09 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:08:10.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_09.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:13.00_UTC" camera=1 focus=10 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:08:13.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_10.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:17.00_UTC" camera=1 focus=11 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:08:17.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_11.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:21.00_UTC" camera=1 focus=12 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:08:21.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_12.bin"
write_flatfield instrument="HMI" flatfield_version=1 series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:25.00_UTC" camera=1 focus=13 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:08:25.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_13.bin" t_stop="2010.04.14_21:05:23.00_UTC"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:28.00_UTC" camera=1 focus=14 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:08:28.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_14.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:32.00_UTC" camera=1 focus=15 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:08:32.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_15.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:36.00_UTC" camera=1 focus=16 t_obs_badpix="2010.04.08_22:07:40.00_UT" t_obs_dark="2010.04.07_14:37:43.00_UTC" t_obs_offpoint="2010.04.08_22:08:36.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_side_16.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:07:38.00_UTC" camera=2 focus=01 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:07:38.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_01.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:07:42.00_UTC" camera=2 focus=02 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:07:42.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_02.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:07:45.00_UTC" camera=2 focus=03 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:07:45.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_03.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:07:49.00_UTC" camera=2 focus=04 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:07:49.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_04.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:07:53.00_UTC" camera=2 focus=05 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:07:53.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_05.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:07:57.00_UTC" camera=2 focus=06 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:07:57.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_06.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:00.00_UTC" camera=2 focus=07 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:08:00.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_07.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:04.00_UTC" camera=2 focus=08 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:08:04.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_08.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:08.00_UTC" camera=2 focus=09 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:08:08.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_09.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:12.00_UTC" camera=2 focus=10 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:08:12.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_10.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:15.00_UTC" camera=2 focus=11 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:08:15.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_11.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:19.00_UTC" camera=2 focus=12 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:08:19.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_12.bin"
write_flatfield instrument="HMI" flatfield_version=1 series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:23.00_UTC" camera=2 focus=13 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:08:23.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_13.bin" t_stop="2010.04.14_21:05:21.00_UTC"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:27.00_UTC" camera=2 focus=14 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:08:27.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_14.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:30.00_UTC" camera=2 focus=15 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:08:30.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_15.bin"
#write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.08_22:08:34.00_UTC" camera=2 focus=16 t_obs_badpix="2010.04.08_22:07:38.00_UT" t_obs_dark="2010.04.07_14:37:41.00_UTC" t_obs_offpoint="2010.04.08_22:08:34.00_UTC" file_flatfield="/tmp20/richard/hmi/flat_front_16.bin"



#//04/.14.2010

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/scr21/richard/hmi/flatfield_side_4062215_inorbit.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.14_21:05:23.00_UTC" fsn_list_offpoint="3883712" focus=13 nx=4096 ny=4096

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/scr21/richard/hmi/flatfield_front_4062215_inorbit.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.14_21:05:21.00_UTC" fsn_list_offpoint="3883713" focus=13 nx=4096 ny=4096

write_dark instrument="HMI" file_dark="/scr21/richard/hmi/dark_side_3878345.bin" series_dark="hmi.dark" camera=1 t_obs="2010.04.14_15:03:13.00_UTC" fsn_list_dark="3878346,3878366,3878386,3878406,3878426" nx=4096 ny=4096

write_dark instrument="HMI" file_dark="/scr21/richard/hmi/dark_front_3878345.bin" series_dark="hmi.dark" camera=2 t_obs='2010.04.14_15:03:11.00_UTC' fsn_list_dark="3878345,3878365,3878385,3878405,3878425" nx=4096 ny=4096

write_badpix instrument="HMI" file_badpix="/scr21/richard/hmi/badpix_front_4062215" series_badpix="hmi.bad_pixel_list" camera=2 t_obs="2010.04.14_21:05:21.00_UTC" nbad=28

write_badpix instrument="HMI" file_badpix="/scr21/richard/hmi/badpix_side_4062215" series_badpix="hmi.bad_pixel_list" camera=1 t_obs="2010.04.14_21:05:23.00_UTC" nbad=43

write_flatfield  flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.14_21:05:23.00_UTC" camera=1 focus=13 t_obs_badpix="2010.04.14_21:05:23.00_UTC" t_obs_dark="2010.04.14_15:03:13.00_UTC" t_obs_offpoint="2010.04.14_21:05:23.00_UTC" file_flatfield="/scr21/richard/hmi/flatfield_side_4062215_inorbit.bin" t_stop="2010.04.30_23:34:36.48_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.14_21:05:21.00_UTC" camera=2 focus=13 t_obs_badpix="2010.04.14_21:05:21.00_UTC" t_obs_dark="2010.04.14_15:03:11.00_UTC" t_obs_offpoint="2010.04.14_21:05:21.00_UTC" file_flatfield="/scr21/richard/hmi/flatfield_front_4062215_inorbit.bin" t_stop="2010.04.30_23:34:38.34_UTC"

#//04/29/2010

write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/scr21/richard/hmi/flatfield_side_4555561_inorbit.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.04.29_20:15:08.00_UTC" fsn_list_offpoint="4555820" nx=4096 ny=4096 focus=11

write_offpoint pztflag=0  instrument="HMI" file_offpoint="/scr21/richard/hmi/flatfield_front_4555561_inorbit.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.04.29_20:15:06.00_UTC" fsn_list_offpoint="4555819" nx=4096 ny=4096 focus=11

write_badpix instrument="HMI" file_badpix="/scr21/richard/hmi/badpix_side_4555561.bin" series_badpix="hmi.bad_pixel_list" camera=1 t_obs="2010.04.29_20:15:08_UT" nbad=45

write_badpix instrument="HMI" file_badpix="/scr21/richard/hmi/badpix_front_4555561.bin" series_badpix="hmi.bad_pixel_list" camera=2 t_obs="2010.04.29_20:15:06.00_UTC" nbad=32

write_dark instrument="HMI" file_dark="/scr21/richard/hmi/dark_side_4555561.bin" series_dark="hmi.dark" camera=1 t_obs="2010.04.29_17:28:43.00_UTC" fsn_list_dark="4555562,4555582,4555602,4555622,4555642" nx=4096 ny=4096

write_dark instrument="HMI" file_dark="/scr21/richard/hmi/dark_front_4555561.bin" series_dark="hmi.dark" camera=2 t_obs="2010.04.29_17:28:41.00_UTC" fsn_list_dark="4555561,4555581,4555601,4555621,4555641" nx=4096 ny=4096

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.30_23:34:36.48_UTC" camera=1 focus=11 t_obs_badpix="2010.04.29_20:15:08.00_UTC" t_obs_dark="2010.04.29_17:28:43.00_UTC" t_obs_offpoint="2010.04.29_20:15:08.00_UTC" file_flatfield="/scr21/richard/hmi/flat_side_20100430.bin" t_stop="2010.05.06_17:28:39.00_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.04.30_23:34:38.34_UTC" camera=2 focus=11 t_obs_badpix="2010.04.29_20:15:06.00_UTC" t_obs_dark="2010.04.29_17:28:41.00_UT" t_obs_offpoint="2010.04.29_20:15:06.00_UTC" file_flatfield="/scr21/richard/hmi/flat_front_20100430.bin" t_stop="2010.05.06_17:28:38.00_UTC"

#//05/06

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.05.06_17:28:39.00_UTC" camera=1 focus=11 t_obs_badpix="2010.04.29_20:15:08.00_UTC" t_obs_dark="2010.04.29_17:28:43.00_UTC" t_obs_offpoint="2010.04.29_20:15:08.00_UTC" file_flatfield="/scr21/richard/hmi/flat_side_20100506.bin" t_stop="2010.05.13_17:35:21.00_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.05.06_17:28:38.00_UTC" camera=2 focus=11 t_obs_badpix="2010.04.29_20:15:06.00_UTC" t_obs_dark="2010.04.29_17:28:41.00_UT" t_obs_offpoint="2010.04.29_20:15:06.00_UTC" file_flatfield="/scr21/richard/hmi/flat_front_20100506.bin" t_stop="2010.05.13_17:35:23.00_UTC"

#//05/13

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.05.13_17:35:21.00_UTC" camera=1 focus=11 t_obs_badpix="2010.04.29_20:15:08.00_UTC" t_obs_dark="2010.04.29_17:28:43.00_UTC" t_obs_offpoint="2010.04.29_20:15:08.00_UTC" file_flatfield="/scr21/richard/hmi/flat_side_20100513.bin" t_stop="2010.05.20_17:37:36.47_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.05.13_17:35:23.00_UTC" camera=2 focus=11 t_obs_badpix="2010.04.29_20:15:06.00_UTC" t_obs_dark="2010.04.29_17:28:41.00_UTC" t_obs_offpoint="2010.04.29_20:15:06.00_UTC" file_flatfield="/scr21/richard/hmi/flat_front_20100513.bin" t_stop="2010.05.20_17:37:38.34_UTC"

#/05/20
write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.05.20_17:37:38.34_UTC" camera=2 focus=11 t_obs_badpix="2010.04.29_20:15:40_TAI" t_obs_dark="2010.04.29_17:29:15_TAI" t_obs_offpoint="2010.04.29_20:15:40_TAI"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_05521253.bin" t_stop="2010.05.27_17:33:04.59_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.05.20_17:37:36.47_UTC" camera=1 focus=11 t_obs_badpix="2010.04.29_20:15:42_TAI" t_obs_dark="2010.04.29_17:29:17_TAI" t_obs_offpoint="2010.04.29_20:15:42_TAI"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_05521253.bin" t_stop="2010.05.27_17:33:02.71_UTC"

#//05/27

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.05.27_17:33:04.59_UTC" camera=2 focus=11 t_obs_badpix="2010.04.29_20:15:40_TAI" t_obs_dark="2010.04.29_17:29:15_TAI" t_obs_offpoint="2010.04.29_20:15:40_TAI"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_05843667.bin" t_stop="2010.06.03_17:42:08.34_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.05.27_17:33:02.71_UTC" camera=1 focus=11 t_obs_badpix="2010.04.29_20:15:42_TAI" t_obs_dark="2010.04.29_17:29:17_TAI" t_obs_offpoint="2010.04.29_20:15:42_TAI"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_05843667.bin"  t_stop="2010.06.03_17:42:06.47_UTC"


#//06.03

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.06.03_17:42:08.34_UTC" camera=2 focus=11 t_obs_badpix="2010.04.29_20:15:40_TAI" t_obs_dark="2010.04.29_17:29:15_TAI" t_obs_offpoint="2010.04.29_20:15:40_TAI"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_06166517.bin" t_stop="2010.06.10_18:06:53_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.06.03_17:42:06.47_UTC" camera=1 focus=11 t_obs_badpix="2010.04.29_20:15:42_TAI" t_obs_dark="2010.04.29_17:29:17_TAI" t_obs_offpoint="2010.04.29_20:15:42_TAI"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_06166517.bin" t_stop="2010.06.10_18:06:51_UTC"



#//06.10
write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_front_06489869.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.06.10_18:06:53_UTC" fsn_list_offpoint="4555819" fsn_list_pzt="6488935,6488937,6488939,6488941,6488943,6488945,6488947,6488949,6488951,6488953,6488955,6489871,6489873,6489875,6489877,6489879,6489881,6489883,6489885,6489887,6489889,6489891" focus=11 nx=4096 ny=4096

write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_side_06489869.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.06.10_18:06:51_UTC" fsn_list_offpoint="4555820" fsn_list_pzt="6488934,6488936,6488938,6488940,6488942,6488944,6488946,6488948,6488950,6488952,6488954,6489870,6489872,6489874,6489876,6489878,6489880,6489882,6489884,6489886,6489888,6489890" focus=11 nx=4096 ny=4096

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.06.10_18:06:53_UTC" camera=2 focus=11 t_obs_badpix="2010.04.29_20:15:40_TAI" t_obs_dark="2010.04.29_17:29:15_TAI" t_obs_offpoint="2010.06.10_18:06:53_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_06489869.bin" t_stop="2010.06.16_18:27:04_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.06.10_18:06:51_UTC" camera=1 focus=11 t_obs_badpix="2010.04.29_20:15:42_TAI" t_obs_dark="2010.04.29_17:29:17_TAI" t_obs_offpoint="2010.06.10_18:06:51_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_06489869.bin" t_stop="2010.06.16_18:27:06_UTC"

#//06.16
write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_front_06766996.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.06.16_18:27:04_UTC" fsn_list_offpoint="4555819"  fsn_list_pzt="6766567,6766569,6766571,6766573,6766575,6766577,6766579,6766581,6766583,6766585,6766587,6766999,6767001,6767003,6767005,6767007,6767009,6767011,6767013,6767015,6767017,6767019" nx=4096 ny=4096 focus=11

write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_side_06766996.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.06.16_18:27:06_UTC" fsn_list_offpoint="4555820" fsn_list_pzt="6766566,6766568,6766570,6766572,6766574,6766576,6766578,6766580,6766582,6766584,6766586,6766998,6767000,6767002,6767004,6767006,6767008,6767010,6767012,6767014,6767016,6767018" nx=4096 ny=4096 focus=11

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.06.16_18:27:04_UTC" camera=2 focus=11 t_obs_badpix="2010.04.29_20:15:40_TAI" t_obs_dark="2010.04.29_17:29:15_TAI" t_obs_offpoint="2010.06.16_18:27:04_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_06766996.bin" t_stop="2010.06.23_18:18:08_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.06.16_18:27:06_UTC" camera=1 focus=11 t_obs_badpix="2010.04.29_20:15:42_TAI" t_obs_dark="2010.04.29_17:29:17_TAI" t_obs_offpoint="2010.06.16_18:27:06_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_06766996.bin" t_stop="2010.06.23_18:18:06_UTC"

#//06.23
write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_front_07089269.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.06.23_18:18:08_UTC" fsn_list_offpoint="4555819"  fsn_list_pzt="7088839,7088841,7088843,7088845,7088847,7088849,7088851,7088853,7088855,7088857,7088859,7089271,7089273,7089275,7089277,7089279,7089281,7089283,7089285,7089287,7089289,7089291" focus=11 nx=4096 ny=4096

write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_side_07089269.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.06.23_18:18:06_UTC" fsn_list_offpoint="4555820" fsn_list_pzt="7088838,7088840,7088842,7088844,7088846,7088848,7088850,7088852,7088854,7088856,7088858,7089270,7089272,7089274,7089276,7089278,7089280,7089282,7089284,7089286,7089288,7089290" focus=11 nx=4096 ny=4096

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.06.23_18:18:06_UTC" camera=1 focus=11 t_obs_badpix="2010.04.29_20:15:42_TAI" t_obs_dark="2010.04.29_17:29:17_TAI" t_obs_offpoint="2010.06.23_18:18:06_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_07089269.bin" t_stop="2010.06.30_18:58:36_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.06.23_18:18:08_UTC" camera=2 focus=11 t_obs_badpix="2010.04.29_20:15:40_TAI" t_obs_dark="2010.04.29_17:29:15_TAI" t_obs_offpoint="2010.06.23_18:18:08_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_07089269.bin" t_stop="2010.06.30_18:58:38_UTC"

#//06.30
write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_side_07413125.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.06.30_18:58:36_UTC" fsn_list_offpoint="4555820" fsn_list_pzt="7412406,7412408,7412410,7412412,7412414,7412416,7412418,7412420,7412422,7412424,7412426,7413126,7413128,7413130,7413132,7413134,7413136,7413138,7413140,7413142,7413144,7413146" nx=4096 ny=4096 focus=11

write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_front_07413125.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.06.30_18:58:38_UTC" fsn_list_offpoint="4555819" fsn_list_pzt="7412407,7412409,7412411,7412413,7412415,7412417,7412419,7412421,7412423,7412425,7412427,7413127,7413129,7413131,7413133,7413135,7413137,7413139,7413141,7413143,7413145,7413147" nx=4096 ny=4096 focus=11

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.06.30_18:58:36_UTC" camera=1 focus=11 t_obs_badpix="2010.04.29_20:15:42_TAI" t_obs_dark="2010.04.29_17:29:17_TAI" t_obs_offpoint="2010.06.30_18:58:36_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_07413125.bin" t_stop="2010.07.07_18:58:36_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.06.30_18:58:38_UTC" camera=2 focus=11 t_obs_badpix="2010.04.29_20:15:40_TAI" t_obs_dark="2010.04.29_17:29:15_TAI" t_obs_offpoint="2010.06.30_18:58:38_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_07413125.bin" t_stop="2010.07.07_18:58:38_UTC"

#//07.07
write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_side_07735685.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.07.07_18:58:36_UTC" fsn_list_offpoint="4555820"  fsn_list_pzt="7735254,7735256,7735258,7735260,7735262,7735264,7735266,7735268,7735270,7735272,7735274,7735686,7735688,7735690,7735692,7735694,7735696,7735698,7735700,7735702,7735704,7735706" nx=4096 ny=4096 focus=11


write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_front_07735685.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.07.07_18:58:38_UTC" fsn_list_offpoint="4555819" fsn_list_pzt="7735255,7735257,7735259,7735261,7735263,7735265,7735267,7735269,7735271,7735273,7735275,7735687,7735689,7735691,7735693,7735695,7735697,7735699,7735701,7735703,7735705,7735707" nx=4096 ny=4096 focus=11


write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.07.07_18:58:36_UTC" camera=1 focus=11 t_obs_badpix="2010.04.29_20:15:42_TAI" t_obs_dark="2010.04.29_17:29:17_TAI" t_obs_offpoint="2010.07.07_18:58:36_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_07735685.bin" t_stop="2010.07.13_19:21:06_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.07.07_18:58:38_UTC" camera=2 focus=11 t_obs_badpix="2010.04.29_20:15:40_TAI" t_obs_dark="2010.04.29_17:29:15_TAI" t_obs_offpoint="2010.07.07_18:58:38_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_07735685.bin" t_stop="2010.07.13_19:21:08_UTC"

#//07.13
write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_side_08012885.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.07.13_19:21:06_UTC" fsn_list_offpoint="4555820" fsn_list_pzt="8012382,8012384,8012386,8012388,8012390,8012392,8012394,8012396,8012398,8012400,8012402,8012886,8012888,8012890,8012892,8012894,8012896,8012898,8012900,8012902,8012904,8012906" nx=4096 ny=4096 focus=11


write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_front_08012885.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.07.13_19:21:08_UTC" fsn_list_offpoint="4555819" fsn_list_pzt="8012383,8012385,8012387,8012389,8012391,8012393,8012395,8012397,8012399,8012401,8012403,8012887,8012889,8012891,8012893,8012895,8012897,8012899,8012901,8012903,8012905,8012907" nx=4096 ny=4096 focus=11

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.07.13_19:21:06_UTC" camera=1 focus=11 t_obs_badpix="2010.04.29_20:15:42_TAI" t_obs_dark="2010.04.29_17:29:17_TAI" t_obs_offpoint="2010.07.13_19:21:06_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_08012885.bin" t_stop="2010.07.21_20:01:36_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.07.13_19:21:08_UTC" camera=2 focus=11 t_obs_badpix="2010.04.29_20:15:40_TAI" t_obs_dark="2010.04.29_17:29:15_TAI" t_obs_offpoint="2010.07.13_19:21:08_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_08012885.bin" t_stop="2010.07.21_20:01:38_UTC"

#///
write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/scr21/richard/hmi/flat_side_11_08097524.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.07.15_22:13:27_TAI" fsn_list_offpoint="8097524" nx=4096 ny=4096 focus=11
write_offpoint  pztflag=0 instrument="HMI" file_offpoint="/scr21/richard/hmi/flat_front_11_08097523.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.07.15_22:13:25_TAI" fsn_list_offpoint="8097523" nx=4096 ny=4096 focus=11

write_dark instrument="HMI" fsn_list_dark="8096854" series_dark="hmi.dark" camera=1 t_obs="2010.07.15_22:19:32_TAI" file_dark="/scr21/richard/hmi/dark_inorbit_side_8108138.bin" nx=4096 ny=4096
write_dark instrument="HMI" fsn_list_dark="8096853" series_dark="hmi.dark" camera=2 t_obs="2010.07.15_22:19:30_TAI" file_dark="/scr21/richard/hmi/dark_inorbit_front_8108138.bin" nx=4096 ny=4096

write_badpix instrument="HMI" file_badpix="/scr21/richard/hmi/badpix_side_4555561.bin" series_badpix="hmi.bad_pixel_list" camera=1 t_obs="2010.07.15_22:13:27_TAI" nbad=45

write_badpix instrument="HMI" file_badpix="/scr21/richard/hmi/badpix_front_4555561.bin" series_badpix="hmi.bad_pixel_list" camera=2 t_obs="2010.07.15_22:13:25_TAI" nbad=32


#//07.21
write_offpoint instrument="HMI"\
file_offpoint="/scr21/richard/hmi/pzt_flat_side_08380193.bin"\
series_offpoint="hmi.offpoint_flatfield" camera=1\
t_obs="2010.07.21_20:01:36_UTC"\
fsn_list_offpoint="8096924"\
fsn_list_pzt="8379762,8379764,8379766,8379768,8379770,8379772,8379774,8379776,8379778,8379780,8379782,8380194,8380196,8380198,8380200,8380202,8380204,8380206,8380208,8380210,8380212,8380214" nx=4096 ny=4096 focus=11

write_offpoint instrument="HMI"\
 file_offpoint="/scr21/richard/hmi/pzt_flat_front_08380193.bin"\
 series_offpoint="hmi.offpoint_flatfield" camera=2\
 t_obs="2010.07.21_20:01:38_UTC"\
 fsn_list_offpoint="8096923"\
fsn_list_pzt="8379763,8379765,8379767,8379769,8379771,8379773,8379775,8379777,8379779,8379781,8379783,8380195,8380197,8380199,8380201,8380203,8380205,8380207,8380209,8380211,8380213,8380215" nx=4096 ny=4096 focus=11

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.07.21_20:01:36_UTC" camera=1 focus=11 t_obs_badpix="2010.07.15_22:13:27_TAI" t_obs_dark="2010.07.15_22:19:32_TAI" t_obs_offpoint="2010.07.21_20:01:36_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_08380193.bin" t_stop="2010.07.28_19:48:06_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.07.21_20:01:38_UTC" camera=2 focus=11 t_obs_badpix="2010.07.15_22:13:25_TAI" t_obs_dark="2010.07.15_22:19:30_TAI" t_obs_offpoint="2010.07.21_20:01:38_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_08380193.bin" t_stop="2010.07.28_19:48:08_UTC"

#//07.28
write_offpoint instrument="HMI"\
 file_offpoint="/scr21/richard/hmi/pzt_flat_side_08702321.bin"\
 series_offpoint="hmi.offpoint_flatfield" camera=1\
 t_obs="2010.07.28_19:48:06_UTC"\
 fsn_list_offpoint="8096924"\
fsn_list_pzt="8701242,8701244,8701246,8701248,8701250,8701252,8701254,8701256,8701258,8701260,8701262,8702322,8702324,8702326,8702328,8702330,8702332,8702334,8702336,8702338,8702340,8702342" nx=4096 ny=4096 focus=11

write_offpoint instrument="HMI"\
 file_offpoint="/scr21/richard/hmi/pzt_flat_front_08702321.bin"\
 series_offpoint="hmi.offpoint_flatfield" camera=2\
 t_obs="2010.07.28_19:48:08_UTC"\
 fsn_list_offpoint="8096923"\
fsn_list_pzt="8701243,8701245,8701247,8701249,8701251,8701253,8701255,8701257,8701259,8701261,8701263,8702323,8702325,8702327,8702329,8702331,8702333,8702335,8702337,8702339,8702341,8702343" nx=4096 ny=4096 focus=11

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.07.28_19:48:06_UTC" camera=1 focus=11 t_obs_badpix="2010.07.15_22:13:27_TAI" t_obs_dark="2010.07.15_22:19:32_TAI" t_obs_offpoint="2010.07.28_19:48:06_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_08702321.bin" t_stop="2010.08.04_19:59:21_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.07.28_19:48:08_UTC" camera=2 focus=11 t_obs_badpix="2010.07.15_22:13:25_TAI" t_obs_dark="2010.07.15_22:19:30_TAI" t_obs_offpoint="2010.07.28_19:48:08_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_08702321.bin" t_stop="2010.08.04_19:59:23_UTC"

#//08.04
write_offpoint instrument="HMI"\
 file_offpoint="/scr21/richard/hmi/pzt_flat_side_09025241.bin"\
 series_offpoint="hmi.offpoint_flatfield" camera=1\
 t_obs="2010.08.04_19:59:21_UTC"\
 fsn_list_offpoint="8096924"\
fsn_list_pzt="9024450,9024452,9024454,9024456,9024458,9024460,9024462,9024464,9024466,9024468,9024470,9025242,9025244,9025246,9025248,9025250,9025252,9025254,9025256,9025258,9025260,9025262" nx=4096 ny=4096 focus=11

write_offpoint instrument="HMI"\
 file_offpoint="/scr21/richard/hmi/pzt_flat_front_09025241.bin"\
 series_offpoint="hmi.offpoint_flatfield" camera=2\
 t_obs="2010.08.04_19:59:23_UTC"\
 fsn_list_offpoint="8096923"\
fsn_list_pzt="9024451,9024453,9024455,9024457,9024459,9024461,9024463,9024465,9024467,9024469,9024471,9025243,9025245,9025247,9025249,9025251,9025253,9025255,9025257,9025259,9025261,9025263" nx=4096 ny=4096 focus=11

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.08.04_19:59:21_UTC" camera=1 focus=11 t_obs_badpix="2010.07.15_22:13:27_TAI" t_obs_dark="2010.07.15_22:19:32_TAI" t_obs_offpoint="2010.08.04_19:59:21_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_09025241.bin" t_stop="2010.08.11_18:56:21_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.08.04_19:59:23_UTC" camera=2 focus=11 t_obs_badpix="2010.07.15_22:13:25_TAI" t_obs_dark="2010.07.15_22:19:30_TAI" t_obs_offpoint="2010.08.04_19:59:23_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_09025241.bin" t_stop="2010.08.11_18:56:23_UTC"

#//08.11
write_offpoint instrument="HMI"\
 file_offpoint="/scr21/richard/hmi/pzt_flat_side_09345785.bin"\
 series_offpoint="hmi.offpoint_flatfield" camera=1\
 t_obs="2010.08.11_18:56:21_UTC"\
 fsn_list_offpoint="8096924"\
fsn_list_pzt="9345210,9345212,9345214,9345216,9345218,9345220,9345222,9345224,9345226,9345228,9345230,9345786,9345788,9345790,9345792,9345794,9345796,9345798,9345800,9345802,9345804,9345806" nx=4096 ny=4096 focus=11

write_offpoint instrument="HMI"\
 file_offpoint="/scr21/richard/hmi/pzt_flat_front_09345785.bin"\
 series_offpoint="hmi.offpoint_flatfield" camera=2\
 t_obs="2010.08.11_18:56:23_UTC"\
 fsn_list_offpoint="8096923"\
fsn_list_pzt="9345211,9345213,9345215,9345217,9345219,9345221,9345223,9345225,9345227,9345229,9345231,9345787,9345789,9345791,9345793,9345795,9345797,9345799,9345801,9345803,9345805,9345807" nx=4096 ny=4096 focus=11

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.08.11_18:56:21_UTC" camera=1 focus=11 t_obs_badpix="2010.07.15_22:13:27_TAI" t_obs_dark="2010.07.15_22:19:32_TAI" t_obs_offpoint="2010.08.11_18:56:21_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_09345785.bin" t_stop="2010.08.18_20:08:21_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.08.11_18:56:23_UTC" camera=2 focus=11 t_obs_badpix="2010.07.15_22:13:25_TAI" t_obs_dark="2010.07.15_22:19:30_TAI" t_obs_offpoint="2010.08.11_18:56:23_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_09345785.bin" t_stop="2010.08.18_20:08:23_UTC"

#//08.18
write_offpoint instrument="HMI"\
file_offpoint="/scr21/richard/hmi/pzt_flat_side_09670469.bin"\
series_offpoint="hmi.offpoint_flatfield" camera=1\
t_obs="2010.08.18_20:08:21_UTC" fsn_list_offpoint="8096924"\
fsn_list_pzt="9669606,9669608,9669610,9669612,9669614,9669616,9669618,9669620,9669622,9669624,9669626,9670470,9670472,9670474,9670476,9670478,9670480,9670482,9670484,9670486,9670488,9670490" nx=4096 ny=4096 focus=11

write_offpoint instrument="HMI"\
 file_offpoint="/scr21/richard/hmi/pzt_flat_front_09670469.bin"\
 series_offpoint="hmi.offpoint_flatfield" camera=2\
 t_obs="2010.08.18_20:08:23_UTC"\
 fsn_list_offpoint="8096923"\
fsn_list_pzt="9669607,9669609,9669611,9669613,9669615,9669617,9669619,9669621,9669623,9669625,9669627,9670471,9670473,9670475,9670477,9670479,9670481,9670483,9670485,9670487,9670489,9670491" nx=4096 ny=4096 focus=11

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.08.18_20:08:21_UTC" camera=1 focus=11 t_obs_badpix="2010.07.15_22:13:27_TAI" t_obs_dark="2010.07.15_22:19:32_TAI" t_obs_offpoint="2010.08.18_20:08:21_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_09670469.bin" t_stop="2010.08.25_17:39:51_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.08.18_20:08:23_UTC" camera=2 focus=11 t_obs_badpix="2010.07.15_22:13:25_TAI" t_obs_dark="2010.07.15_22:19:30_TAI" t_obs_offpoint="2010.08.18_20:08:23_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_09670469.bin" t_stop="2010.08.25_17:39:53_UTC"


#//08.25
write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_side_09987737.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.08.25_17:39:51_UTC" fsn_list_offpoint="8096924" fsn_list_pzt="9987162,9987164,9987166,9987168,9987170,9987172,9987174,9987176,9987178,9987180,9987182,9987738,9987740,9987742,9987744,9987746,9987748,9987750,9987752,9987754,9987756,9987758" nx=4096 ny=4096 focus=11

write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_front_09987737.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.08.25_17:39:53_UTC" fsn_list_offpoint="8096923" fsn_list_pzt="9987163,9987165,9987167,9987169,9987171,9987173,9987175,9987177,9987179,9987181,9987183,9987739,9987741,9987743,9987745,9987747,9987749,9987751,9987753,9987755,9987757,9987759" nx=4096 ny=4096 focus=11

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.08.25_17:39:51_UTC" camera=1 focus=11 t_obs_badpix="2010.07.15_22:13:27_TAI" t_obs_dark="2010.07.15_22:19:32_TAI" t_obs_offpoint="2010.08.25_17:39:51_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_09987737.bin" t_stop="2010.09.01_20:10:36_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.08.25_17:39:53_UTC" camera=2 focus=11 t_obs_badpix="2010.07.15_22:13:25_TAI" t_obs_dark="2010.07.15_22:19:30_TAI" t_obs_offpoint="2010.08.25_17:39:53_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_09987737.bin" t_stop="2010.09.01_20:10:38_UTC"

#//09.01
write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_side_10315121.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.09.01_20:10:36_UTC" fsn_list_offpoint="8096924" fsn_list_pzt="10314186,10314188,10314190,10314192,10314194,10314196,10314198,10314200,10314202,10314204,10314206,10315122,10315124,10315126,10315128,10315130,10315132,10315134,10315136,10315138,10315140,10315142" nx=4096 ny=4096 focus=11

write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_front_10315121.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.09.01_20:10:38_UTC" fsn_list_offpoint="8096923" fsn_list_pzt="10314187,10314189,10314191,10314193,10314195,10314197,10314199,10314201,10314203,10314205,10314207,10315123,10315125,10315127,10315129,10315131,10315133,10315135,10315137,10315139,10315141,10315143" nx=4096 ny=4096 focus=11

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.09.01_20:10:36_UTC" camera=1 focus=11 t_obs_badpix="2010.07.15_22:13:27_TAI" t_obs_dark="2010.07.15_22:19:32_TAI" t_obs_offpoint="2010.09.01_20:10:36_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_10315121.bin" t_stop="2010.09.08_19:18:51_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.09.01_20:10:38_UTC" camera=2 focus=11 t_obs_badpix="2010.07.15_22:13:25_TAI" t_obs_dark="2010.07.15_22:19:30_TAI" t_obs_offpoint="2010.09.01_20:10:38_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_10315121.bin" t_stop="2010.09.08_19:18:53_UTC"

#//09.08
write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_side_10636025.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.09.08_19:18:51_UTC" fsn_list_offpoint="8096924" fsn_list_pzt="10635378,10635380,10635382,10635384,10635386,10635388,10635390,10635392,10635394,10635396,10635398,10636026,10636028,10636030,10636032,10636034,10636036,10636038,10636040,10636042,10636044,10636046" nx=4096 ny=4096 focus=11

write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_front_10636025.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.09.08_19:18:53_UTC" fsn_list_offpoint="8096923" fsn_list_pzt="10635379,10635381,10635383,10635385,10635387,10635389,10635391,10635393,10635395,10635397,10635399,10636027,10636029,10636031,10636033,10636035,10636037,10636039,10636041,10636043,10636045,10636047" nx=4096 ny=4096 focus=11

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.09.08_19:18:51_UTC" camera=1 focus=11 t_obs_badpix="2010.07.15_22:13:27_TAI" t_obs_dark="2010.07.15_22:19:32_TAI" t_obs_offpoint="2010.09.08_19:18:51_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_10636025.bin" t_stop="2010.09.15_19:07:36_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.09.08_19:18:53_UTC" camera=2 focus=11 t_obs_badpix="2010.07.15_22:13:25_TAI" t_obs_dark="2010.07.15_22:19:30_TAI" t_obs_offpoint="2010.09.08_19:18:53_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_10636025.bin" t_stop="2010.09.15_19:07:38_UTC"


#//09.15
write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_side_10958225.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.09.15_19:07:36_UTC" fsn_list_offpoint="8096924" fsn_list_pzt="10957794,10957796,10957798,10957800,10957802,10957804,10957806,10957808,10957810,10957812,10957814,10958226,10958228,10958230,10958232,10958234,10958236,10958238,10958240,10958242,10958244,10958246" nx=4096 ny=4096 focus=11

write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_front_10958225.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.09.15_19:07:38_UTC" fsn_list_offpoint="8096923" fsn_list_pzt="10957795,10957797,10957799,10957801,10957803,10957805,10957807,10957809,10957811,10957813,10957815,10958227,10958229,10958231,10958233,10958235,10958237,10958239,10958241,10958243,10958245,10958247" nx=4096 ny=4096 focus=11

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.09.15_19:07:36_UTC" camera=1 focus=11 t_obs_badpix="2010.07.15_22:13:27_TAI" t_obs_dark="2010.07.15_22:19:32_TAI" t_obs_offpoint="2010.09.15_19:07:36_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_10958225.bin" t_stop="2010.09.22_18:51:51_UTC"

write_flatfield flatfield_version=1 instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.09.15_19:07:38_UTC" camera=2 focus=11 t_obs_badpix="2010.07.15_22:13:25_TAI" t_obs_dark="2010.07.15_22:19:30_TAI" t_obs_offpoint="2010.09.15_19:07:38_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_10958225.bin" t_stop="2010.09.22_18:51:51_UTC"

#//09.22
write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_side_11280281.bin" series_offpoint="hmi.offpoint_flatfield" camera=1 t_obs="2010.09.22_18:51:51_UTC" fsn_list_offpoint="8096924" fsn_list_pzt="11279850,11279852,11279854,11279856,11279858,11279860,11279862,11279864,11279866,11279868,11279870,11280282,11280284,11280286,11280288,11280290,11280292,11280294,11280296,11280298,11280300,11280302" nx=4096 ny=4096 focus=11

write_offpoint instrument="HMI" file_offpoint="/scr21/richard/hmi/pzt_flat_front_11280281.bin" series_offpoint="hmi.offpoint_flatfield" camera=2 t_obs="2010.09.22_18:51:53_UTC" fsn_list_offpoint="8096923" fsn_list_pzt="11279851,11279853,11279855,11279857,11279859,11279861,11279863,11279865,11279867,11279869,11279871,11280283,11280285,11280287,11280289,11280291,11280293,11280295,11280297,11280299,11280301,11280303" nx=4096 ny=4096 focus=11

write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.09.22_18:51:51_UTC" camera=1 focus=11 t_obs_badpix="2010.07.15_22:13:27_TAI" t_obs_dark="2010.07.15_22:19:32_TAI" t_obs_offpoint="2010.09.22_18:51:51_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_side_11280281.bin"

write_flatfield instrument="HMI" series_offpoint="hmi.offpoint_flatfield" series_badpix="hmi.bad_pixel_list" series_dark="hmi.dark" series_flatfield="hmi.flatfield" t_start="2010.09.22_18:51:53_UTC" camera=2 focus=11 t_obs_badpix="2010.07.15_22:13:25_TAI" t_obs_dark="2010.07.15_22:19:30_TAI" t_obs_offpoint="2010.09.22_18:51:53_UTC"  file_flatfield="/scr21/richard/hmi/pzt_flat_front_11280281.bin"


