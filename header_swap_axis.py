
'''
Swap the axes in a header, without losing keys to WCS
'''

from astropy.wcs import WCS


def header_swapaxes(header, ax1, ax2):
    '''
    '''

    mywcs = WCS(header)

    new_hdr = mywcs.swapaxes(ax1, ax2).to_header()

    lost_keys = list(set(header.keys()) - set(new_hdr.keys()))

    for key in lost_keys:
        # Drop COMMENT, HISTORY, "" for now as they are a pain
        if key == "" or key == "COMMENT" or key == "HISTORY":
            continue

        if str(ax1+1) in key:
            new_hdr[key.replace(str(ax1+1), str(ax2+1))] = header[key]
        elif str(ax2+1) in key:
            new_hdr[key.replace(str(ax2+1), str(ax1+1))] = header[key]
        else:
            new_hdr[key] = header[key]

    return new_hdr
