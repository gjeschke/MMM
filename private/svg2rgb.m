function rgb = svg2rgb(svgcolor)
% rgb = svg2rgb(svgcolor)
% returns the rgb values for a SVG color name (normalized)
% an empty variable is returned if the color is undefined

switch lower(svgcolor)
    case 'aliceblue'
        rgb = [240,248,255];
    case 'antiquewhite'
        rgb = [250,235,215];
    case 'aqua'
        rgb = [0,255,255];
    case 'aquamarine'
        rgb = [127,255,212];
    case 'azure'
        rgb = [240,255,255];
    case 'beige'
        rgb = [245,245,220];
    case 'bisque'
        rgb = [255,228,196];
    case 'black'
        rgb = [0,0,0];
    case 'blanchedalmond'
        rgb = [255,235,205];
    case 'blue'
        rgb = [0,0,255];
    case 'blueviolet'
        rgb = [138,43,226];
    case 'brown'
        rgb = [165,42,42];
    case 'burlywood'
        rgb = [222,184,135];
    case 'cadetblue'
        rgb = [95,158,160];
    case 'chartreuse'
        rgb = [127,255,0];
    case 'chocolate'
        rgb = [210,105,30];
    case 'coral'
        rgb = [255,127,80];
    case 'cornflowerblue'
        rgb = [100,149,237];
    case 'cornsilk'
        rgb = [255,248,220];
    case 'crimson'
        rgb = [220,20,60];
    case 'cyan'
        rgb = [0,255,255];
    case 'darkblue'
        rgb = [0,0,139];
    case 'darkcyan'
        rgb = [0,139,139];
    case 'darkgoldenrod'
        rgb = [184,134,11];
    case 'darkgray'
        rgb = [169,169,169];
    case 'darkgreen'
        rgb = [0,100,0];
    case 'darkgrey'
        rgb = [169,169,169];
    case 'darkkhaki'
        rgb = [189,183,107];
    case 'darkmagenta'
        rgb = [139,0,139];
    case 'darkolivegreen'
        rgb = [85,107,47];
    case 'darkorange'
        rgb = [255,140,0];
    case 'darkorchid'
        rgb = [153,50,204];
    case 'darkred'
        rgb = [139,0,0];
    case 'darksalmon'
        rgb = [233,150,122];
    case 'darkseagreen'
        rgb = [143,188,143];
    case 'darkslateblue'
        rgb = [72,61,139];
    case 'darkslategray'
        rgb = [47,79,79];
    case 'darkslategrey'
        rgb = [47,79,79];
    case 'darkturquoise'
        rgb = [0,206,209];
    case 'darkviolet'
        rgb = [148,0,211];
    case 'deeppink'
        rgb = [255,20,147];
    case 'deepskyblue'
        rgb = [0,191,255];
    case 'dimgray'
        rgb = [105,105,105];
    case 'dimgrey'
        rgb = [105,105,105];
    case 'dodgerblue'
        rgb = [30,144,255];
    case 'firebrick'
        rgb = [178,34,34];
    case 'floralwhite'
        rgb = [255,250,240];
    case 'forestgreen'
        rgb = [34,139,34];
    case 'fuchsia'
        rgb = [255,0,255];
    case 'gainsboro'
        rgb = [220,220,220];
    case 'ghostwhite'
        rgb = [248,248,255];
    case 'gold'
        rgb = [255,215,0];
    case 'goldenrod'
        rgb = [218,165,32];
    case 'gray'
        rgb = [128,128,128];
    case 'green'
        rgb = [0,128,0];
    case 'greenyellow'
        rgb = [173,255,47];
    case 'grey'
        rgb = [128,128,128];
    case 'honeydew'
        rgb = [240,255,240];
    case 'hotpink'
        rgb = [255,105,180];
    case 'indianred'
        rgb = [205,92,92];
    case 'indigo'
        rgb = [75,0,130];
    case 'ivory'
        rgb = [255,255,240];
    case 'khaki'
        rgb = [240,230,140];
    case 'lavender'
        rgb = [230,230,250];
    case 'lavenderblush'
        rgb = [255,240,245];
    case 'lawngreen'
        rgb = [124,252,0];
    case 'lemonchiffon'
        rgb = [255,250,205];
    case 'lightblue'
        rgb = [173,216,230];
    case 'lightcoral'
        rgb = [240,128,128];
    case 'lightcyan'
        rgb = [224,255,255];
    case 'lightgoldenrodyellow'
        rgb = [250,250,210];
    case 'lightgray'
        rgb = [211,211,211];
    case 'lightgreen'
        rgb = [144,238,144];
    case 'lightgrey'
        rgb = [211,211,211];
    case 'lightpink'
        rgb = [255,182,193];
    case 'lightsalmon'
        rgb = [255,160,122];
    case 'lightseagreen'
        rgb = [32,178,170];
    case 'lightskyblue'
        rgb = [135,206,250];
    case 'lightslategray'
        rgb = [119,136,153];
    case 'lightslategrey'
        rgb = [119,136,153];
    case 'lightsteelblue'
        rgb = [176,196,222];
    case 'lightyellow'
        rgb = [255,255,224];
    case 'lime'
        rgb = [0,255,0];
    case 'limegreen'
        rgb = [50,205,50];
    case 'linen'
        rgb = [250,240,230];
    case 'magenta'
        rgb = [255,0,255];
    case 'maroon'
        rgb = [128,0,0];
    case 'mediumaquamarine'
        rgb = [102,205,170];
    case 'mediumblue'
        rgb = [0,0,205];
    case 'mediumorchid'
        rgb = [186,85,211];
    case 'mediumpurple'
        rgb = [147,112,219];
    case 'mediumseagreen'
        rgb = [60,179,113];
    case 'mediumslateblue'
        rgb = [123,104,238];
    case 'mediumspringgreen'
        rgb = [0,250,154];
    case 'mediumturquoise'
        rgb = [72,209,204];
    case 'mediumvioletred'
        rgb = [199,21,133];
    case 'midnightblue'
        rgb = [25,25,112];
    case 'mintcream'
        rgb = [245,255,250];
    case 'mistyrose'
        rgb = [255,228,225];
    case 'moccasin'
        rgb = [255,228,181];
    case 'navajowhite'
        rgb = [255,222,173];
    case 'navy'
        rgb = [0,0,128];
    case 'oldlace'
        rgb = [253,245,230];
    case 'olive'
        rgb = [128,128,0];
    case 'olivedrab'
        rgb = [107,142,35];
    case 'orange'
        rgb = [255,165,0];
    case 'orangered'
        rgb = [255,69,0];
    case 'orchid'
        rgb = [218,112,214];
    case 'palegoldenrod'
        rgb = [238,232,170];
    case 'palegreen'
        rgb = [152,251,152];
    case 'paleturquoise'
        rgb = [175,238,238];
    case 'palevioletred'
        rgb = [219,112,147];
    case 'papayawhip'
        rgb = [255,239,213];
    case 'peachpuff'
        rgb = [255,218,185];
    case 'peru'
        rgb = [205,133,63];
    case 'pink'
        rgb = [255,192,203];
    case 'plum'
        rgb = [221,160,221];
    case 'powderblue'
        rgb = [176,224,230];
    case 'purple'
        rgb = [128,0,128];
    case 'red'
        rgb = [255,0,0];
    case 'rosybrown'
        rgb = [188,143,143];
    case 'royalblue'
        rgb = [65,105,225];
    case 'saddlebrown'
        rgb = [139,69,19];
    case 'salmon'
        rgb = [250,128,114];
    case 'sandybrown'
        rgb = [244,164,96];
    case 'seagreen'
        rgb = [46,139,87];
    case 'seashell'
        rgb = [255,245,238];
    case 'sienna'
        rgb = [160,82,45];
    case 'silver'
        rgb = [192,192,192];
    case 'skyblue'
        rgb = [135,206,235];
    case 'slateblue'
        rgb = [106,90,205];
    case 'slategray'
        rgb = [112,128,144];
    case 'slategrey'
        rgb = [112,128,144];
    case 'snow'
        rgb = [255,250,250];
    case 'springgreen'
        rgb = [0,255,127];
    case 'steelblue'
        rgb = [70,130,180];
    case 'tan'
        rgb = [210,180,140];
    case 'teal'
        rgb = [0,128,128];
    case 'thistle'
        rgb = [216,191,216];
    case 'tomato'
        rgb = [255,99,71];
    case 'turquoise'
        rgb = [64,224,208];
    case 'violet'
        rgb = [238,130,238];
    case 'wheat'
        rgb = [245,222,179];
    case 'white'
        rgb = [255,255,255];
    case 'whitesmoke'
        rgb = [245,245,245];
    case 'yellow'
        rgb = [255,255,0];
    case 'yellowgreen'
        rgb = [154,205,50];
    otherwise
        rgb = [];
end
rgb = rgb/255;


