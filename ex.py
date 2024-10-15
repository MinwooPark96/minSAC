from src.model import CausalDiagram as CD
from src.msbd import MultiOutcomeSequentialBackdoorCriterion as mSBD
from src.sac import SequentialAdjustmentCriterion as SAC


G_1 = CD(
    vs = {'X1','X2'} | {'Y1','Y2'} | {'Za','Zb','Zc','Zd'},
    directed_edges =(('Zb', 'Y1'),
    ('Zb', 'X2'),
    ('Zb', 'X1'),
    ('Zd', 'Y2'),
    ('Za', 'Zb'),
    ('Za', 'X1'),
    ('Zc', 'Zd'),
    ('Zc', 'X2'),
    ('X2', 'Y2'),
    ('X1', 'Y1'),
    ('Y1', 'X2')),
    bidirected_edges = [
    ('Zb', 'Zc', 'U0')
    ])



if __name__ == '__main__':
    xs,ys,zs = {'X1','X2'},{'Y1','Y2'},{'Za','Zb','Zc','Zd'}
    msbd = mSBD(G_1,xs,ys,zs)
    msbd.set_covariate(1,{'Za','Zb'})
    msbd.set_covariate(2,{'Zc','Zd'})
    
    print('mSBD criterion : ',msbd.criteria())
    for e in msbd.G.edges:
        if e in  msbd.localGraph(1).edges:
            pass
        else :
            print(e, 'is removed')

    print('-------------------')

    for e in msbd.G.edges:
        if e in  msbd.localGraph(2).edges:
            pass
        else :
            print(e, 'is removed')

    print('-------------------')

    # msbd.formula()

    # (mSBD -> SAC)
    sac = SAC(G_1,xs,ys,zs)

    sac.set_covariate(1,{'Za','Zb'})
    sac.set_covariate(2,{'Zc','Zd'})

    for e in sac.G.edges:
        if e in  sac.psbd(1).edges:
            pass
        else :
            print(e, 'is removed')

    print('-------------------')


    for e in sac.G.edges:
        if e in  sac.psbd(2).edges:
            pass
        else :
            print(e, 'is removed')

    print('-------------------')

