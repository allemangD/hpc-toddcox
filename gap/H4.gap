f := FreeGroup( 4 );
g := f / [ f.1^2,
           f.2^2,
           f.3^2,
           f.4^2,
           ( f.1*f.2 )^5,
           ( f.2*f.3 )^3,
           ( f.3*f.4 )^3,
           ( f.1*f.3 )^2,
           ( f.1*f.4 )^2,
           ( f.2*f.4 )^2
         ];

tab := CosetTable( g, Subgroup(g, [ ]));
time;
