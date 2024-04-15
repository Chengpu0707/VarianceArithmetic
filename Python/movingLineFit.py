import enum

import varDbl

class FitType (enum.Enum):
    LOCAL = 1,
    MOVING = 2,
    MOVING_NO_VAR_ADJ = -1,



def movingLineFit(sInput, half:int, fitType:FitType) ->tuple[tuple[tuple[varDbl.VarDbl]]]:
    '''
    For each "sInput" timeseries of length n, line-fit in the window with half size ("half"*2 + 1)
    according to "fitType":
     *) LOCAL: not progressively, as a testing standard
     *) MOVING_ADJUST: progressive moving window with adjusting the variance
     *) MOVING_NO_ADJUST: progressive moving window without adjusting the variance,
    '''
    full = 2 * half + 1
    denom0 = 1 /full
    denom1 = 3 /half /(half + 1) /full

    def value(val):
        return val.value() if type(val) == varDbl.VarDbl else val

    match fitType:
        case FitType.LOCAL:
            return tuple([(sum(sInput[j - half: j + half + 1])*denom0,
                           sum([k * sInput[j + k] for k in range(-half, half + 1)])*denom1)
                          for j in range(half, len(sInput) - half)])
        case FitType.MOVING_NO_VAR_ADJ:
            c0 = sum(sInput[1 : full])
            c1 = sum([k * sInput[half + k] for k in range(-half, half + 1)])
            sFit = [((sInput[0] + c0)*denom0, c1*denom1)]
            for j in range(full, len(sInput)):
                c1 += half*(sInput[j - full] + sInput[j]) - c0
                c0 += sInput[j] - sInput[j - full + 1]
                sFit.append(((sInput[j - full + 1] + c0)*denom0, c1*denom1))
            return tuple(sFit)
        
        case FitType.MOVING:
            denom0sq = denom0**2
            denom1sq = denom1**2
            c0 = sum([value(sInput[k]) for k in range(1, full)])
            c1 = sum([k * value(sInput[half + k])for k in range(-half, half + 1)])
            v0 = sum([sInput[k].variance() for k in range(full) if type(sInput[k]) == varDbl.VarDbl])
            v1 = sum([k**2 * sInput[half + k].variance()
                      for k in range(-half, half + 1) if type(sInput[half + k]) == varDbl.VarDbl])
            sFit = [(varDbl.VarDbl((value(sInput[0]) + c0)*denom0, v0 *denom0sq, True),
                     varDbl.VarDbl(c1 *denom1, v1 *denom1sq, True))]
            for j in range(full, len(sInput)):
                c1 += half*(value(sInput[j - full]) + value(sInput[j])) - c0
                v1 = sum([k**2 * sInput[j - half + k].variance()
                          for k in range(-half, half + 1) if type(sInput[j - half + k]) == varDbl.VarDbl])
                
                c0 += value(sInput[j]) - value(sInput[j - full + 1])
                if (type(sInput[j]) == varDbl.VarDbl):
                    v0 += sInput[j].variance()
                if (type(sInput[j - full + 1]) == varDbl.VarDbl):
                    v0 -= sInput[j - full + 1].variance()

                sFit.append((varDbl.VarDbl((value(sInput[j - full + 1]) + c0)*denom0, v0 *denom0sq, True), 
                             varDbl.VarDbl(c1 *denom1, v1 *denom1sq, True)))
            return tuple(sFit)
        
        case _:
            raise NotImplemented()
    
