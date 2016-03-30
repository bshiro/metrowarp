var ImgWarper = ImgWarper || {};

ImgWarper.Warper = function(orig, gridDim, optGridSize, optAlpha) {
	this.alpha = optAlpha || 1;
  this.gridSize = optGridSize || 20;
  this.width = gridDim[0];
  this.height = gridDim[1];
  //this.bilinearInterpolation = new ImgWarper.BilinearInterpolation(this.width, this.height, canvas);

	this.cellRows = 0;
	this.cellCols = 0;
	
  this.grid = [];
	this.transformedGrid = [];
  for (var i = 0; i <= this.width ; i += this.gridSize) {
		this.cellRows = 0;
    for (var j = 0; j <= this.height ; j += this.gridSize) {
			
      a = new ImgWarper.Point(i,j);
      b = new ImgWarper.Point(i + this.gridSize, j);
      c = new ImgWarper.Point(i + this.gridSize, j + this.gridSize);
      d = new ImgWarper.Point(i, j + this.gridSize);
      this.grid.push([a, b, c, d]);
			this.cellRows++;
    }
		this.cellCols++;
  }
	
	this.oriPoints = new Array();
  this.dstPoints = new Array();	
	for (var i=0; i<orig.length; i++)
	{
		var q = new ImgWarper.Point(orig[i].x, orig[i].y);
		this.oriPoints.push(q);
		this.dstPoints.push(q);
	}


}

ImgWarper.Warper.prototype.warp = function(fromPoints, toPoints) {
  var deformation = new ImgWarper.AffineDeformation(toPoints, fromPoints, this.alpha);
  for (var i = 0; i < this.grid.length; ++i) {
    this.transformedGrid[i] = [
        deformation.pointMover(this.grid[i][0]),
        deformation.pointMover(this.grid[i][1]),
        deformation.pointMover(this.grid[i][2]),
        deformation.pointMover(this.grid[i][3])];
  }

}

ImgWarper.AffineDeformation = function(fromPoints, toPoints, alpha) {
  this.w = null;
  this.pRelative = null;
  this.qRelative = null;
  this.A = null;
  if (fromPoints.length != toPoints.length) {
    console.error('Points are not of same length.'); 
    return;
  }
  this.n = fromPoints.length;  
  this.fromPoints = fromPoints;
  this.toPoints = toPoints;
  this.alpha = alpha;
};

ImgWarper.AffineDeformation.prototype.pointMover = function (point){

  if (null == this.pRelative || this.pRelative.length < this.n) {
    this.pRelative = new Array(this.n); 
  }
  if (null == this.qRelative || this.qRelative.length < this.n) {
    this.qRelative = new Array(this.n); 
  }
  if (null == this.w || this.w.length < this.n) {
    this.w = new Array(this.n);
  }
  if (null == this.A || this.A.length < this.n) {
    this.A = new Array(this.n); 
  }

  for (var i = 0; i < this.n; ++i) {
    var t = this.fromPoints[i].subtract(point);
    this.w[i] = Math.pow(t.x * t.x + t.y * t.y, -this.alpha);
  }
  var pAverage = ImgWarper.Point.weightedAverage(this.fromPoints, this.w);
  var qAverage = ImgWarper.Point.weightedAverage(this.toPoints, this.w);

  for (var i = 0; i < this.n; ++i) {
    this.pRelative[i] = this.fromPoints[i].subtract(pAverage);
    this.qRelative[i] = this.toPoints[i].subtract(qAverage);
  }
  var B = new ImgWarper.Matrix22(0, 0, 0, 0);

  for (var i = 0; i < this.n; ++i) {
    B.addM(this.pRelative[i].wXtX(this.w[i]));
  }

  B = B.inverse();
  for (var j = 0; j < this.n; ++j) {
    this.A[j] = point.subtract(pAverage).multiply(B)
      .dotP(this.pRelative[j]) * this.w[j];
  }

  var r = qAverage; //r is an point 
  for (var j = 0; j < this.n; ++j) {
    r = r.add(this.qRelative[j].multiply_d(this.A[j]));
  }
  //return r;
	return point.multiply_d(2).subtract(r);
};
