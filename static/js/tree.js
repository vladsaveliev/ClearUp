// Based on http://bl.ocks.org/mbostock/4063570

function buildTree(svgElemId, treeWidth, treeHeight, data) {
    var treeSampleSpace = 8,
        sampleSeqSpace = 15;

    var svg = d3.select(svgElemId);
    svg.attr("width", treeWidth);
    svg.attr("height", treeHeight);
    var g = svg.append("g").attr("transform", "translate(40,0)");

    var tree = d3.cluster()
        .size([treeHeight, treeWidth]);

    var stratify = d3.stratify()
        .parentId(function(d) { return d.id.substring(0, d.id.lastIndexOf(".")); });

    var root = stratify(data)
        .sort(function(a, b) { return (a.height - b.height) || a.id.localeCompare(b.id); });

    tree(root);

    var link = g.selectAll(".link")
        .data(root.descendants().slice(1))
      .enter().append("path")
        .attr("class", "link")
        .attr("d", function(d) {
          return "M" + d.y + "," + d.x
               + "C" + (d.parent.y + 0) + "," + d.x
               + " " + (d.parent.y + 10) + "," + d.parent.x
               + " " + d.parent.y + "," + d.parent.x;
        });

    g.selectAll(".node")
        .data(root.descendants())
      .enter().append("g")
        .attr("class", function(d) { return "node" + (d.children ? " node--internal" : " node--leaf"); })
        .attr("transform", function(d) { return "translate(" + d.y + "," + d.x + ")"; });

    var leafsTexts = g.selectAll(".node--leaf")
      .append("a")
        .attr("xlink:href", function(d) { return d.data.sample_id; })
      .append("text")
        .attr("dy", 3)
        .attr("x", treeSampleSpace)
        .style("text-anchor", "start")
        .style("fill", function(d) { return d.data.color; })
        .text(function(d) { return d.id.substring(d.id.lastIndexOf(".") + 1); });

    var labelSizes = [];
    leafsTexts.each(function() { labelSizes.push(this.getComputedTextLength()); });
    var maxTextWidth = d3.max(labelSizes);

    var seqWidth = g.selectAll(".node--leaf")
      .append("text")
        .attr("dy", 3)
        .attr("x", function() { return maxTextWidth + sampleSeqSpace; })
        .attr("class", "seq")
        .each(function(d, i) {
          d3.select(this)
            .selectAll(".nuc")
              .data(d.data.seq)
            .enter().append("tspan")
              .attr("class", function(d) { return d[1]; })
              .text(function(d) { return d[0]; })
            .append("tspan")
              .attr("width", (i > 0 && i % 2 == 0) ? 3 : 0 );
        })
        .node()
        .getComputedTextLength();

    svg.attr("width",
        treeWidth +
        treeSampleSpace +
        maxTextWidth +
        sampleSeqSpace +
        seqWidth +
        50);

    console.log(seqWidth);
}
