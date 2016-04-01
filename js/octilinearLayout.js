/**
 * Created by pio on 3/30/2016.
 */
var OctolinearLayout = OctolinearLayout || {};

OctolinearLayout.conjgrad = function(A, b, x) {
    var Ax = numeric.dotMV(A,x);    //A*x
    var r = numeric.sub(b,Ax);      //r=b-Ax
    var p = r.slice(0);     //shallow copy
    var rsold = numeric.dotVV(r,r);//rsold=r'*r
    for (var i=0; i< b.length; i++)
    {
        var Ap = numeric.dotMV(A,p);                 //Ap=A*p;
        var alpha = rsold / numeric.dotVV(p,Ap);     //alpha=rsold/(p'*Ap);
        numeric.addeq(x, numeric.mul(alpha, p));    //x=x+alpha*p;
        numeric.subeq(r, numeric.mul(alpha, Ap));   //r=r-alpha*Ap;
        var rsnew = numeric.dotVV(r,r);             //rsnew=r'*r
        if (Math.sqrt(rsnew)<1e-10)
            break;
        p = numeric.add(r, numeric.mul(rsnew/rsold, p));    //p=r+(rsnew/rsold)*p;
        rsold = rsnew;
    }
    return x;
}

OctolinearLayout.getAngle = function(v1, v2){
    //return Math.acos( (v1[0] * v2[0] + v1[1] * v2[1]) / ( Math.sqrt(v1[0]*v1[0] + v1[1]*v1[1]) * Math.sqrt(v2[0]*v2[0] + v2[1]*v2[1]) ) );
    return Math.atan2((v1[0] * v2[1] - v1[1] * v2[0]), v1[0] * v2[0] + v1[1] * v2[1]);
}

OctolinearLayout.modLengthObjective = function(bx, by, x, y, vertices, edges, weight, meanEdgeLength){
    var nEdges = edges.length;
    for (var eid=0; eid<nEdges; eid++) {
        var i = edges[eid][0];
        var j = edges[eid][1];

        var v_ij = numeric.sub(vertices[i],vertices[j]);
        var v_hat_ij = numeric.sub([x[i],y[i]],[x[j],y[j]]);
        var theta = -this.getAngle(v_hat_ij,v_ij);

        //var R = [[Math.cos(edgeTheta[eid]), -Math.sin(edgeTheta[eid])] ,[Math.sin(edgeTheta[eid]), Math.cos(edgeTheta[eid])]];
        var R = [[Math.cos(theta), -Math.sin(theta)] ,[Math.sin(theta), Math.cos(theta)]];
        var s_ij = meanEdgeLength / numeric.norm2(v_ij);
        var b = numeric.dotMV(numeric.mul(s_ij, R),v_ij);

        bx[eid] = weight*b[0];
        by[eid] = weight*b[1];
    }
}

OctolinearLayout.computeEnergy = function(A, x, b){
    return numeric.norm2Squared(numeric.sub(numeric.dotMV(A,x),b));
}

OctolinearLayout.appendBoundaryObjective = function(A, b, x, xMax) {
    var cols = x.length;
    for (var i=0; i< cols; i++)
    {
        if (x[i]<0) {
            var row = new Array(cols).fill(0);
            row[i] = 1;
            A.push(row);
            b.push(0);
            //console.log("append");
        } else if (x[i]>xMax) {
            var row = new Array(cols).fill(0);
            row[i] = 1;
            A.push(row);
            b.push(xMax);
            //console.log("append");
        }
    }
}

OctolinearLayout.curvilinearLayout = function(vertices, edges, distortedStations) {
    //test conjgrad
    //var x = this.conjgrad([[4,1],[1,3]], [1,2], [2,1]);
    //console.log(x);
    //console.log(vertices);
    //console.log(edges);
    //console.log(this.getAngle(numeric.sub([5,2],[1,1]),numeric.sub([10,10],[2,2])));
   // console.log(this.getAngle(numeric.sub([10,10],[2,2]),numeric.sub([5,2],[1,1])));

    var weights = {"maxangle":1, "length":5, "guide":0.05}

    var nStations = vertices.length;
    var nEdges = edges.length;

    var neighbors = [];
    var x = [];
    var y = [];
    for (var i = 0; i < nStations; i++) {
        neighbors.push([]);
        x.push(vertices[i][0]);
        y.push(vertices[i][1]);
    }
    var meanEdgeLength = 0;
    for (var i = 0; i < edges.length; i++) {
        var v_i = edges[i][0];
        var v_j = edges[i][1];
        neighbors[v_i].push(v_j);
        neighbors[v_j].push(v_i);

        meanEdgeLength += numeric.norm2(numeric.sub(vertices[v_i], vertices[v_j]));
    }
    meanEdgeLength /= nEdges;

    var bx = [];
    var A = [];
    var by = [];
    //var Ay = [];



    var edgeTheta = new Array(nEdges).fill(0);
    for (var eid = 0; eid < nEdges; eid++) {
        var i = edges[eid][0];
        var j = edges[eid][1];
        var R = [[Math.cos(edgeTheta[eid]), -Math.sin(edgeTheta[eid])], [Math.sin(edgeTheta[eid]), Math.cos(edgeTheta[eid])]];
        var v_ij = numeric.sub(vertices[i], vertices[j]);
        var s_ij = meanEdgeLength / numeric.norm2(v_ij);
        var b = numeric.dotMV(numeric.mul(s_ij, R), v_ij);

        var row = new Array(nStations).fill(0);
        row[i] = weights["length"]*1;
        row[j] = weights["length"]*-1;

        A.push(row);
        bx.push(weights["length"]*b[0]);
        by.push(weights["length"]*b[1]);
    }


    for (var v_i = 0; v_i < nStations; v_i++) {
        var f = neighbors[v_i].length;
        if (f<2)
            continue;
        var theta = 2 * Math.PI / f;
        var tanCoef = Math.tan((Math.PI - theta) * 0.5);

        for (var j = 0; j < f; j++) {
            var row = new Array(nStations).fill(0);
            var k = (j + 1) % f;
            var v_j = neighbors[v_i][j];
            var v_k = neighbors[v_i][k];
            //console.log(v_i+","+v_j+","+v_k);
            row[v_i] = 1;
            row[v_j] = -1 + 0.5 + tanCoef * 0.5;
            row[v_k] = -0.5 - tanCoef * 0.5;

            A.push(row);
            //console.log(row);
            bx.push(0);
            by.push(0);
        }
    }

    for (var v_i = 0; v_i < nStations; v_i++) {
        var row = new Array(nStations).fill(0);
        row[v_i] = weights["guide"]*1;
        A.push(row);
        bx.push(weights["guide"]*x[v_i]);
        by.push(weights["guide"]*y[v_i]);
    }

    var k = 0;
    var maxIter = 50;

    var Ax = A.slice(0);
    var Ay = A.slice(0);
    var xMax = Math.max(...x);
    var yMax = Math.max(...y);

    var energy = this.computeEnergy(Ax,x,bx)+this.computeEnergy(Ay,y,by);
    var prevEnergy = energy;
    var energydiff = (prevEnergy - energy)/prevEnergy;
    do {
        var Atx = numeric.transpose(Ax);
        var Aty = numeric.transpose(Ay);
        x = this.conjgrad(numeric.dotMMbig(Atx, Ax), numeric.dotMV(Atx, bx), x);
        y = this.conjgrad(numeric.dotMMbig(Aty, Ay), numeric.dotMV(Aty, by), y);

        energy = this.computeEnergy(Ax,x,bx)+this.computeEnergy(Ay,y,by);
        //console.log(energy);

        this.modLengthObjective(bx, by, x, y, vertices, edges, weights["length"], meanEdgeLength);
        this.appendBoundaryObjective(Ax,bx,x,xMax);
        this.appendBoundaryObjective(Ay,by,y,yMax);
        k++;

        energydiff = (prevEnergy - energy)/prevEnergy;
    }while ( k<maxIter && energydiff>1e-4 );

    for (var i=0; i<nStations; i++)
    {
        distortedStations[i].x = x[i];
        distortedStations[i].y = y[i];
    }
};

OctolinearLayout.octLayout = function(vertices, edges, distortedStations) {

    var nStations = vertices.length;
    var nEdges = edges.length;

    var meanEdgeLength = 0;
    for (var i = 0; i < edges.length; i++) {
        var v_i = edges[i][0];
        var v_j = edges[i][1];
        meanEdgeLength += numeric.norm2(numeric.sub(vertices[v_i], vertices[v_j]));
    }
    meanEdgeLength /= nEdges;

    var weights = {"octilinear":10}

    var x = [];
    var y = [];
    for (var i = 0; i < nStations; i++) {
        x.push(distortedStations[i].x);
        y.push(distortedStations[i].y);
    }

    var bx = [];
    var A = [];
    var by = [];

    var pi4th = Math.PI/4;
    for (var eid = 0; eid < nEdges; eid++) {
        var i = edges[eid][0];
        var j = edges[eid][1];

        var v_ij = numeric.sub([x[i],y[i]], [x[j],y[j]]);
        var horVec = [meanEdgeLength,0];
        //var theta = Math.round(this.getAngle(horVec,v_ij) / (pi4th)) * pi4th;
        var v_ij_angle = OctolinearLayout.getAngle(horVec,v_ij);
        var theta = (Math.round(v_ij_angle / pi4th) * pi4th)-v_ij_angle;

        var R = [[Math.cos(theta), -Math.sin(theta)] ,[Math.sin(theta), Math.cos(theta)]];
        var b = numeric.dotMV(R, v_ij);
        var row = new Array(nStations).fill(0);
        row[i] = weights["octilinear"]*1;
        row[j] = weights["octilinear"]*-1;

        A.push(row);
        bx.push(weights["octilinear"]*b[0]);
        by.push(weights["octilinear"]*b[1]);
    }
    for (var v_i = 0; v_i < nStations; v_i++) {
        var row = new Array(nStations).fill(0);
        row[v_i] = 1;
        A.push(row);
        bx.push(x[v_i]);
        by.push(y[v_i]);
    }

   // return;
    var xMax = Math.max(...x);
    var yMax = Math.max(...y);
    var k = 0;
    var maxIter = 100;
    var Ax = A.slice(0);
    var Ay = A.slice(0);
    var energy = this.computeEnergy(Ax,x,bx)+this.computeEnergy(Ay,y,by);
    var prevEnergy = energy;
    var energydiff = (prevEnergy - energy)/prevEnergy;
    do {
        var Atx = numeric.transpose(Ax);
        var Aty = numeric.transpose(Ay);
        x = this.conjgrad(numeric.dotMMbig(Atx, Ax), numeric.dotMV(Atx, bx), x);
        y = this.conjgrad(numeric.dotMMbig(Aty, Ay), numeric.dotMV(Aty, by), y);

        energy = this.computeEnergy(Ax,x,bx)+this.computeEnergy(Ay,y,by);
        this.appendBoundaryObjective(Ax,bx,x,xMax);
        this.appendBoundaryObjective(Ay,by,y,yMax);
        k++;

        energydiff = (prevEnergy - energy)/prevEnergy;
    }while ( k<maxIter && energydiff>1e-4 );

    for (var i=0; i<nStations; i++)
    {
        distortedStations[i].x = x[i];
        distortedStations[i].y = y[i];
    }
};

OctolinearLayout.lineSegmentIntersection = function(a,b, c,d) {
    var u = numeric.sub(b,a);
    var v = numeric.sub(d,c);
    var w = numeric.sub(a,c);
    var den = -v[0] * u[1] + u[0] * v[1];
    var s = (-u[1] * w[0] + u[0] * w[1]) / den;
    var t = (v[0] * w[1] - v[1] * w[0]) / den;
    return (s >= 0 && s <= 1 && t >= 0 && t <= 1);
};

OctolinearLayout.intersectionCheck = function(x,y,edges){
    //naive method
    intersections = [];
    for (var i=0; i<edges.length; i++) {
        var e_i0 = [x[edges[i][0]], y[edges[i][0]]];
        var e_i1 = [x[edges[i][1]], y[edges[i][1]]];
        for (var i=0; i<edges.length; i++) {
            if (i==j)
                continue;
            var e_j0 = [x[edges[j][0]], y[edges[j][0]]];
            var e_j1 = [x[edges[j][1]], y[edges[j][1]]];

            if (this.lineSegmentIntersection(e_i0,e_i1,e_j0, e_j1))
                intersections.push([i,j]);

        }
    }

};