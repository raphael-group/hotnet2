d3.json("hotnet_viz_data.json", function(error, samples) {

  if(error) {
    console.log('error: ');
    console.log(error);
    return;
  }

  var animation_speed = 300;
  var coloring = {"BLCA": "#A6CEE3",
                    "BRCA": "#1F78B4",
                    "COADREAD":"#B2DF8A",
                    "GBM":"#33A02C",
                    "HNSC":"#FB9A99",
                    "KIRC":"#E31A1C",
                    "LAML":"#FDBF6F",
                    "LUAD":"#FF7F00",
                    "LUSC":"#CAB2D6",
                    "OV":"#6A3D9A",
                    "READ":"#FFFF99",
                    "UCEC":"#B15928"},
      highlightColor = "#f1c40f",
      selectedColor = "#E74C3C",
      textColorStrongest = "#2C3E50",
      textColorStrong = "#34495E",
      textColorLight = "#BDC3C7",
      textColorLightest = "#ECF0F1"
      blockColorStrongest = "#2C3E50"
      blockColorStrong = "#7F8C8D"
      blockColorMedium = "#95A5A6",
      blockColorLight = "#BDC3C7",
      blockColorLightest = "#ECF0F1",
      heatRed = "#E74C3C",
      heatBlue = "#3498DB";

  var networkDimensions = { "width": 500, "height" :520},
      oncoprintDimensions = { "width": 800, "height" :800},
      infoboxDimensions = { "width": 500, "height" :400},
      annotationDimensions = { "width": 800, "height": 500};

  var gene_display_list;

  var menu = d3.select("body")
    .append("div");

  samples.forEach(function(sample){
    createButtons(menu.append("div").attr("class", "dashboard"), sample)
  })

  function createButtons(DOMcontainer, data){

    DOMcontainer.append("ul")
      .selectAll("li")
      .data(data.nodes).enter()
      .append("li")
      .style("background-color", blockColorMedium)
      .text(function(d){ return d.gene; })
      .style("color", textColorLightest);

    var canvas = DOMcontainer.append("div")
      .attr("class", "unselected")
      .attr("id", function(d, i){
        return "number" + i;
      });

    DOMcontainer.append("p")
      .attr("class", "toggleButton")
      .on("click", function(d, i){ 
        toggleButton(DOMcontainer, canvas, data); 
      })
      .text("...");
  }

  function toggleButton(DOMcontainer, canvas, data_set){
    if (canvas.classed("unselected")){
      canvas.transition().duration(animation_speed/2)
        .style("height", function(){
          return oncoprintDimensions["height"] + 1000+ "px"
        })
        .style("width", function(){
          return oncoprintDimensions["width"] + networkDimensions["width"] + 200 +  "px"
        });

      setTimeout(function(){
        generateDashboard(DOMcontainer, canvas, data_set);
        }, animation_speed);
    }
    else{
      canvas.selectAll("div").remove();
      canvas.attr("class", "unselected");
      canvas.transition().duration(animation_speed)
        .style("height", "0px");
      }
  }

  function generateDashboard(DOMcontainer, canvas, data_set){

    var display_dict = {
      "gene": geneList = d3.keys(countGenes(data_set.samples.slice(), "gene")),
      "cancer": cancerList = d3.keys(countCategory(data_set.samples.slice(), "cancer")),
      "mutation": mutationList = d3.keys(countGenes(data_set.samples.slice(), "CNA"))
    }

    var min = d3.min(data_set.nodes, function(d){ return d.heat;}),
        max = d3.max(data_set.nodes, function(d){ return d.heat;}),
        heat_color_scale = d3.scale.linear()
          .domain([min, max])
          .range([heatBlue, heatRed]);


    heatColorScale = function heatColorScale(heatValue){
      return heat_color_scale(heatValue);
    };

    toggleCategory = function toggleCategory(datum, category){
      var updated_list = updateList(datum, category);
      gene_display_list = updated_list;

      network.updateGenes();
      oncoprint.updateGenes();
      loliplot.updateGenes();

      DOMcontainer.selectAll("li")
        .filter(function(d){ return gene_display_list.indexOf(d.gene) > -1 })
        .transition().duration(animation_speed/2)
        .style("color", textColorLightest)

      DOMcontainer.selectAll("li")
        .filter(function(d){ return gene_display_list.indexOf(d.gene) <= -1 })
        .transition().duration(animation_speed/2)
        .style("color", blockColorMedium)
      
    }

    highlightCategory = function highlightCategory(datum, category){
//      highlightGenes(network.node.selectAll(".node"),function(d){ return d.gene; }, datum, "class", "highlighted");
//      highlightGenes( network.selectAll(".node"),function(d){ return d.gene; }, datum, "stroke", highlightColor);

      network.highlightGenes(datum);
      oncoprint.highlightGenes(datum);
      loliplot.highlightGenes(datum);
      DOMcontainer.selectAll("li")
        .filter(function(d){ return d.gene == datum })
        .transition().duration(animation_speed/2)
        .style("background-color", blockColorStrongest)
    }

    unhighlightCategory = function unhighlightCategory(){
      network.unhighlightGenes();
      oncoprint.unhighlightGenes();
      loliplot.unhighlightGenes();
      DOMcontainer.selectAll("li")
        .transition().duration(animation_speed/2)
        .style("background-color", blockColorLight)
    }

    updateList = function updateList(datum, category){
      
      var list = display_dict[category];
      if (list.indexOf(datum) > -1){
        var removed = list.splice(list.indexOf(datum), 1);
      }
      else{
        list.push(datum);
      }
      return list;
    }

    gene_display_list = display_dict["gene"].slice();

    var leftColumn = canvas.append("div")
                          .attr("id", "leftColumn")
                          .attr("class", "column")
                          .style("width", networkDimensions["width"] + "px"),
        rightColumn = canvas.append("div")
                          .attr("id", "rightColumn")
                          .attr("class", "column")
                          .style("width", oncoprintDimensions["width"] + "px");

    var network = new generateNetwork(leftColumn.append("div").attr("class", "network"), data_set, this);
    var oncoprint = new generateOncoprint(rightColumn.append("div").attr("class", "oncoprint"), data_set, this);
    var infobox = new generateInfoBox(leftColumn.append("div").attr("class", "infobox"), data_set, this);
    var loliplot = new generateLoliplot(rightColumn.append("div").attr("class", "loliplot"), data_set, this);
    DOMcontainer.selectAll("li")
      .transition().duration(animation_speed/2)
      .style("background-color", blockColorLight)
    canvas.attr("class", "selected");

    makeInteractive(this, DOMcontainer.selectAll("li"), function(d){ return d.gene; }, "gene");

    oncoprint.updateGenes(display_dict["gene"].slice());

  }

  function generateInfoBox(DOMcontainer, data_set, dashboard){
    var width = infoboxDimensions["width"],
        height = infoboxDimensions["height"];
    DOMcontainer
      .attr("width", width + "px")
      .attr("height", height + "px")
      .append("p")
      .text(" Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.");
  }

  function generateNetwork(DOMcontainer, data_set, dashboard){

    var width = networkDimensions["width"],
        height = networkDimensions["height"];

    var numNodes = data_set.nodes.length,
        numLinks = data_set.links.filter(function(d){ return d.present }).length,
        linkStrength = (height)*numLinks/(numNodes);

    var force = d3.layout.force()
        .charge(-3000)
        .linkDistance(200)
        .size([width, width]);

    DOMcontainer.append("h1").text("HOTNET NETWORK");

    DOMcontainer
      .attr("width", width)
      .attr("height", height);

    var svg = DOMcontainer
      .append("div")
      .append("svg") 
      .attr("width", width)
      .attr("height", width);

    var legend = DOMcontainer
      .append("div")
      .append("svg") 
      .attr("width", width)
      .attr("height", height - width);

    var legendBlockWidth = (height-width),
        numBlocks = Math.floor(width/legendBlockWidth),
        min = d3.min(data_set.nodes, function(d){ return d.heat;}),
        max = d3.max(data_set.nodes, function(d){ return d.heat;}),
        interval = (max - min) / (numBlocks);

    legend.selectAll("rect")
      .data(function(){
        var temp_array = [];
        for (var i = 0; i < numBlocks; i++){
          temp_array.push(i*interval);
        }
        return temp_array;
      }).enter()
      .append("rect")
      .attr("width", legendBlockWidth)
      .attr("height", legendBlockWidth)
      .attr("x", function(d, i){ return i*legendBlockWidth})
      .attr("y", 0)
      .attr("fill", function(d, i){return heatColorScale(min + i*interval)});

    legend.selectAll("text").data([min, max]).enter().append("text")
      .attr("x", function(d, i){ return  i*width})
      .attr("y", (height-width) - 2)
      .attr("text-anchor", function(d, i){
        return (i == 0)? "start":"end" ;
      })
      .text(function(d){ 
        return String((d*1000).toFixed(2)) + " e-3"
      })
      .attr("font-size", 14)
      .attr("fill", textColorLightest);
    force.nodes(data_set.nodes)
      .links(data_set.links)

    var link = svg.selectAll(".link")
      .data(data_set.links)
      .enter()
      .append("line")
      .on("mouseover", function(d){
        d3.select(this)
          .transition().duration(animation_speed/2)
          .style("stroke", highlightColor)
        dashboard.highlightCategory(d.source.gene, "gene");
        dashboard.highlightCategory(d.target.gene, "gene");
      })
      .on("mouseout", function(d){
        d3.select(this)
          .transition().duration(animation_speed/2)
          .style("stroke", blockColorLight)
        dashboard.unhighlightCategory();
      })
      .attr("class", "link");

    link.filter(function(d){ return !d.present })
      .style("stroke-width", 0);

    var node = svg.selectAll(".node")
      .data(data_set.nodes)
      .enter().append("g")
      .attr("class", "node");

    makeInteractive(dashboard, node, function(d){ return d.gene; }, "gene");

    var count_dict = countGenes(data_set.samples, "gene");

    node.append("circle")
//      .attr("class", "circle")
      .attr("r", function(d, i){
        return (count_dict[d.gene]/count_dict["total"]*100)
      });

    node.style("fill", function(d) { 
        return dashboard.heatColorScale(d.heat); 
      })
      .append("text")
      .attr("x", function(d, i){
        return (count_dict[d.gene]/count_dict["total"]*100) + 12
      })
      .attr("dy", ".35em")
      .text(function(d) { return d.gene; });

    force.on("tick", function() {
      link.attr("x1", function(d) { return d.source.x; })
          .attr("y1", function(d) { return d.source.y; })
          .attr("x2", function(d) { return d.target.x; })
          .attr("y2", function(d) { return d.target.y; });
      node.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });
    });

    force.start();
    for (var i = 0; i < 100; ++i){ 
      force.tick();
    }

    this.highlightGenes = highlightGenes;
    this.unhighlightGenes = unhighlightGenes;
    this.updateGenes = updateGenes;
    this.heatColorScale = heatColorScale;

    function updateGenes(){
      var hideNodes = node.filter(function(d){
          return gene_display_list.indexOf(d.gene) <= -1;
        });

      hideNodes.style("fill-opacity", 0.25);

      hideNodes.selectAll("circle")  
        .transition().duration(animation_speed)
        .attr("r", 5);

      var keepNodes = node.filter(function(d){
          return gene_display_list.indexOf(d.gene) > -1;
        });

      keepNodes.style("fill-opacity", 1);

      keepNodes.selectAll("circle")
        .transition().duration(animation_speed)
        .attr("r", function(d, i){
          return count_dict[d.gene]/count_dict["total"]*100;
        });

      var hideLinks = link.filter(function(d){
        return gene_display_list.indexOf(d.source.gene) <= -1 || gene_display_list.indexOf(d.target.gene) <= -1
      })
      .filter(function(d){
        return d.present;
      })
      .transition().duration(animation_speed)
      .style("stroke-opacity", 0.25)

      var keepLinks = link.filter(function(d){
        return gene_display_list.indexOf(d.source.gene) > -1 && gene_display_list.indexOf(d.target.gene) > -1
      })
      .filter(function(d){
        return d.present;
      })
      .transition().duration(animation_speed)
      .style("stroke-opacity", 1)
    }

    function highlightGenes(gene){

      node.filter(function(d){
        return d.gene == gene;
      })
      .selectAll("circle")
//      .transition().duration(animation_speed)
      .attr("r", function(d, i){
        var tempRad = (gene_display_list.indexOf(d.gene) > -1) ? count_dict[d.gene]/count_dict["total"]*100 : 0;
        return (tempRad > 15) ? tempRad : 15
      });

      var highlight = node.filter(function(d, i){
          return d.gene == gene && gene_display_list.indexOf(gene) > -1;
        });

      highlight.selectAll("circle")
        .transition().duration(animation_speed/2)
        .style("stroke", highlightColor)
        .style("stroke-width", 5)
        .style("stroke-opacity", 1);
    }

    function unhighlightGenes(){
      
      node.selectAll("circle").style("stroke-width", 0);

      node
        .selectAll("circle")
//        .transition().duration(animation_speed)
        .attr("r", function(d, i){
          return (gene_display_list.indexOf(d.gene) > -1) ? count_dict[d.gene]/count_dict["total"]*100 : 5;

        })
        .style("stroke-opacity", 1);
    }
  }

  function generateLoliplot(DOMcontainer, data_set, dashboard){

    DOMcontainer.append("h1").text("GENE ANNOTATIONS");

    var numGenes = data_set.annotations.length,
        boxMargin = 5,
        radius = 5,
        boxWidth = annotationDimensions["width"],
        boxHeight = 120,
        graphWidth = boxWidth,
        graphHeight = (numGenes + 1)*boxMargin + numGenes*(boxHeight+50),
        width = graphWidth,
        height = graphHeight,
        resolution = Math.floor(graphWidth/(radius*2))
        bar_y = 20;  

    var symboling = {"Nonsense_Mutation": 0, "Frame_Shift_Del": 1, "Frame_Shift_Ins": 1, "Missense_Mutation": 2, "Splice_Site": 3, "In_Frame_Del": 4, "In_Frame_Ins": 4}

    var sample = sortViaCategory(data_set.annotations.slice(), "transcripts");

    sample.forEach(function(s){
      if (s.transcripts.length > 0){
        loliButtons(DOMcontainer, s, "h2")
      }
    });


    this.highlightGenes = function highlightGenes(gene){
      DOMcontainer.selectAll(".lolibox")
      .filter(function(d){
        return d == gene;
      })
      .selectAll(".background")
      .style("fill", blockColorStrongest)
    }

    this.unhighlightGenes = function unhighlightGenes(){
      
      DOMcontainer.selectAll(".lolibox")
      .filter(function(d){
        return gene_display_list.indexOf(d) > -1;
      })
      .selectAll(".background")
      .style("fill", blockColorLightest)

      DOMcontainer.selectAll(".lolibox")
      .filter(function(d){
        return gene_display_list.indexOf(d) <= -1;
      })
      .selectAll(".background")
      .style("fill", blockColorLight)
    }

    this.updateGenes = function updateGenes(){

      var hide = DOMcontainer.selectAll(".lolibox")
        .filter(function(d){
          return gene_display_list.indexOf(d) <= -1;
        })

      hide.selectAll(".symbols").style("fill-opacity", 0.25)
          .style("stroke-opacity", 0.25);

      hide.selectAll(".background").style("fill", blockColorLight)

      var keep = DOMcontainer.selectAll(".lolibox")
        .filter(function(d){
          return gene_display_list.indexOf(d) > -1;
        });

      keep.selectAll(".symbols")
          .style("fill-opacity", 1)
          .style("stroke-opacity", 1);

    }

    function loliButtons(selection, sample, DOMElement){
      var div = selection.append("div");

      div.selectAll(".null")
        .data(function(){
          if (sample.transcripts){
            return [sample.gene]
          }
          else{
            return sample
          }
        }).enter()
        .append(DOMElement)
        .text(function(d){ return (sample.transcripts)? d : d.name; })
        .attr("class", "unselected")
        .on("click", function(d){

          if (d3.select(this).classed("unselected")){
            // console.log(sample)
            if (sample.transcripts){
              loliButtons(div,sample.transcripts, "h3");
            }
            else{
              drawLoliplot(div, d, d.name)
            }
            d3.select(this).attr("class", "selected");
          }
          else{
            if (sample.transcripts){
              div.selectAll("svg").remove();
              div.selectAll("h3").remove();
            }
            else{
              div.selectAll("#" + d.name).remove();
            }            
            d3.select(this).attr("class", "unselected");
          }
        });

//      makeInteractive(dashboard, div.selectAll("h2"), function(d){ return d; }, "gene");
    }

    function drawLoliplot(selection, sample, id){
      // var min = d3.min(sample.mutations, function(gene){ return gene.locus;}),
      //     max = d3.max(sample.mutations, function(gene){ return gene.locus;});
      var min = 0,
          max = sample.length;

      var x = d3.scale.linear()
            .domain([min, max])
            .range([0, graphWidth]);

      var xAxis = d3.svg.axis()
            .scale(x)
            .orient("bottom")
            .ticks(5)
            .tickSize(0)
            .tickPadding(1.25);

      var zoom = d3.behavior.zoom()
            .x(x)
            .scaleExtent([1, 20])
            .on("zoom", function(){
              draw(svg);
            });

      var svg = selection
        .selectAll(".null").data([sample.gene]).enter()
        .append("svg")
        .attr("class", "lolibox")
        .attr("id", id)
        .style("width", boxWidth)
        .style("height", boxHeight)
        .call(zoom);

      var background = svg.append("rect")
          .attr("width", boxWidth)
          .attr("height", boxHeight)
          .attr("class", "background")
          .style("fill", blockColorLightest);

      var menu_bar = svg.append("rect")
          .attr("class", "menu_bar")
          .attr("y", boxHeight - bar_y*2)
          .attr("x", x(min))
          .attr("width", x(max)-x(min) - 2*boxMargin)
          .attr("height", bar_y - boxMargin)
          .style("fill", blockColorLight);

      svg.append("g")
        .attr("class", "xaxis")
        .attr("transform", "translate(0," + (boxHeight - bar_y + boxMargin) + ")")
        .style("font-size", "12px")
        .style("fill", blockColorLight)
        .call(xAxis);

      var domains = svg.selectAll(".domains")
        .data(sample.domains.slice()).enter()
        .append("rect")
        .attr("width", function(d, i){
          return x(d.end) - x(d.start);
        })
        .attr("height", bar_y + boxMargin)
        .attr("y", boxHeight - bar_y*2 - boxMargin)
        .attr("class", "domains")
        .style("fill-opacity", .5)
        .style("fill", blockColorMedium);

      var gene_name = sample.gene;
      var data_set = sample["mutations"].slice();
      
      var circle = svg.selectAll(".symbols")
        .data(data_set).enter()
        .append("path")
        .attr("class", "symbols")
        .attr("d", d3.svg.symbol()
          .type(function(d, i){
            return d3.svg.symbolTypes[symboling[d.mutation]];
          })
          .size(radius*radius)
        )
        .style("stroke", function(d, i){
          return coloring[d.cancer];
        })
        .style("fill", function(d, i){
          return coloring[d.cancer];
        })
        .style("stroke-width", 2);

      draw(svg);

      makeInteractive(dashboard, svg, function(d){ return d; }, "gene");

      function draw(selection){
        var cur_min = d3.min(x.domain()),
            cur_max = d3.max(x.domain()),
            cur_res = Math.round((cur_max - cur_min)/resolution);
        cur_res = (cur_res) ? cur_res : 1;

        var index_dict = {};
        for (var i = Math.floor(cur_min/cur_res) - 5; i < Math.ceil(cur_max/cur_res) + 5; i++){
          index_dict[i] = 0;
        }
        var selected_circles = selection.selectAll(".symbols")
          .attr("transform", function(d, i){
            var cur_index = (Math.round(d.locus/cur_res));
            index_dict[cur_index] ++;
            // console.log(x(cur_index*cur_res) + ", " + (bar_y - index_dict[cur_index] * radius*2) );
            return "translate(" + x(cur_index*cur_res) + ", " + (boxHeight - bar_y*2 - index_dict[cur_index] * radius*2 - boxMargin) + ")";
          })
          .attr("stroke-opacity", 1)
          .attr("fill-opacity", 1);

        selection.selectAll(".symbols").filter(function(d, i){ 
            return !(cur_min < d.locus && cur_max > d.locus); 
          })
        .attr("stroke-opacity", 0)
        .attr("fill-opacity", 0);

//        selection.selectAll(".line")
//          .attr("d",line);

       selection.select(".xaxis").call(xAxis);
       selection.select(".menu_bar")
         .attr("x", x(min) + boxMargin)
         .attr("width", x(max)-x(min) - 2*boxMargin)

        selection.selectAll(".domains")
          .attr("x", function(d, i){
            return x(d.start);
          })
          .attr("width", function(d, i){
            return x(d.end) - x(d.start);
          });
      }
    }

  }

  function generateOncoprint(DOMcontainer, data_set, dashboard){
    DOMcontainer.append("h1").text("ONCOPRINTS");
    var width = oncoprintDimensions["width"],
        height = oncoprintDimensions["height"],
        labelHeight = 0,
        labelWidth = 120
        boxMargin = 5,
        graphWidth = width,
        graphHeight = height;

    var numSamples = data_set.samples.length,
        numGenes = data_set.nodes.length;   // We're cheating by using the node length
    var boxWidth,
        boxHeight;

    var x = d3.scale.linear()
      .domain([0, numSamples])
      .range([(labelWidth + boxMargin), (width  - 2*boxMargin)]);

    // Iterates through all samples to organize
    // all genes by size, and returns a dictionary
    // associating each gene to its index
    var samples = data_set.samples.slice(),
        gene_index;

    updateDimensions();

    function updateDimensions(){
      numGenes = gene_display_list.length;
      var cur_min = d3.min(x.domain()),
          cur_max = d3.max(x.domain()),
          tempBoxHeight = (graphHeight - labelHeight - ((numGenes+1)*boxMargin) )/(numGenes);
      boxWidth = (graphWidth - labelWidth - 2*boxMargin)/(cur_max - cur_min);
      boxHeight = 40;
      height = (numGenes)*(boxHeight+boxMargin) + boxMargin;
      gene_index = createGeneOrder(samples);
    }
    function updateGraph(){
      DOMcontainer.transition()
      .duration(animation_speed/2)
      .style("width", width).style("height", height);
      svg.transition()
      .duration(animation_speed/2)
      .style("width", width)
      .style("height", height)
    }
    var zoom = d3.behavior.zoom()
      .x(x)
      .scaleExtent([1, Math.round(20*numSamples/graphWidth)])
      .on("zoom", function(){
        renderOncoprint();
      });

    var svg = DOMcontainer
      .append("svg")
      .call(zoom);

    updateGraph()

    samples.forEach(function(d){
      var cur_genes = d.genes,
          i = Infinity,
          g;
      cur_genes.forEach(function(d){
        if (gene_index[d.gene] < i){
          i = gene_index[d.gene];
          g = d.gene;
        }
      });
      d.coocurrence = cur_genes.length;
      d.primaryGene = g;
    });

    // Takes in sample list, list of gene order, and depth or recursion
    var semisorted_samples = sortViaCategory(samples, "cancer")
    var sorted_samples = sortViaCoocurrence(semisorted_samples, 
      d3.keys(gene_index).sort(function(a, b){ 
        return gene_index[a] - gene_index[b]; 
      }), 1);

    var labels = svg.selectAll(".geneLabels")
      .data(d3.keys(gene_index)
        .sort(function(a, b){return gene_index[a] - gene_index[b] }))
      .enter()
      .append("svg:g")
      .attr("class", "geneLabels");

    labels.append("rect")
      .attr("class", "geneLabelsRow")
      .attr("fill", blockColorLight);

    // Takes in the dashboard, selection of things to make
    // interactive, and the category of what we're interacting
    // with (i.e. cancer, genes, mutation type)
    makeInteractive(dashboard, labels, function(d){ return d; }, "gene");

    var data_box = svg.append("svg:g");
    
    var boxes = data_box.selectAll(".oncosample")
      .data(sorted_samples)
      .enter().append("svg:g")
      .attr("class", "oncosample")
      .attr("fill", function(d){
        return coloring[d.cancer]
      })
      .attr("id", function(d){
        return d.name;
      })
      .on("mouseover", function(d){
        highlightSample(d3.select(this), true);
      })
      .on("mouseout", function(d){
        unhilightSample(d3.select(this), true);
      });

    boxes.append("text")
      .attr("text-anchor", "end")
      .text(function(d){
        var nameLen = d.name.length;
        return d.name.slice(nameLen - 4, nameLen);
      })
      .attr("fill", blockColorMedium);



    // boxes.selectAll(".box")
    //   .enter(function(d){ return d.genes }).enter()
    //     .append("rect")
    //     .attr("class", "box");

    sorted_samples.forEach(function(d, i){
      data_box.selectAll("#" + d.name)
        .selectAll(".box")
        .data(d.genes).enter()
        .append("rect")
        .attr("class", "box");
    })


    var count_dict = countGenes(samples, "gene");

    labels.append("text")
      .attr("class", "geneLabelsText")
      .attr("font-size", 14)
      .attr("fill", textColorLightest)
      .text(function(d){
        return d + " (" + count_dict[d] + ")";
      });

    boxes.selectAll(".inactivating")
    .data(function(d){ 
      return d.genes; })
    .enter()
    .append("rect")
    .filter(function(d){ return d.inactivating})
    .attr("class", "inactivating")
    .style("fill", blockColorStrongest)
    .attr("width", boxWidth)
    .attr("height", boxHeight/4)


    this.renderOncoprint = renderOncoprint;
    this.highlightGenes = highlightGenes;
    this.unhighlightGenes = unhighlightGenes;
    this.updateGenes = updateGenes;

    function renderOncoprint(){
      updateDimensions();
      updateGraph();
      renderBoxes(boxes, gene_index);
      renderLabels(labels, gene_index);
    }

    function updateGenes(){

      renderOncoOpacities(svg.selectAll(".geneLabels").selectAll(".geneLabelsRow"), 1, 1)

      data_box.selectAll(".oncosample")
        .filter(function(d){
          return !checkValidGene(d.genes, gene_display_list) 
        })
        .style("fill-opacity", 0)
        .selectAll(".box")
        .transition().duration(animation_speed/2)
        .style("fill-opacity", 0);

      boxes = data_box.selectAll(".oncosample")
        .filter(function(d){
          return checkValidGene(d.genes, gene_display_list) });

      renderOncoprint();
    }

    function highlightGenes(gene){
      var highlight = boxes
        .filter(function(d, i){
          var valid = false;
          d.genes.forEach(function(g){
            if (g.gene == gene){
              valid = true;
            }
          })
          return valid;
        });
      highlightSample(highlight, false);

      labels.filter(function(d, i){
          return d == gene && gene_display_list.indexOf(gene) > -1;
        })
        .selectAll("rect")
        .transition().duration(animation_speed/2)
        .style("fill", blockColorStrongest)
    }

    function unhighlightGenes(){
      unhilightSample(data_box.selectAll(".oncosample"), false)

     labels.selectAll(".geneLabelsRow")
      .transition().duration(animation_speed/2)
      .style("fill", "#BDC3C7")
    }

    function highlightSample(selection, boolText){
      selection
      //.selectAll(".box")
        .filter(function(d){ return gene_display_list.indexOf(d.gene) > -1})
        .transition().duration(animation_speed/2)
        .style("fill", highlightColor);

      if (boolText){
        selection.selectAll("text")
          .transition().duration(animation_speed/2)
          .style("font-size", 10)
          .style("fill", "#2C3E50")        
      }
    }

    function unhilightSample(selection, boolText){
      selection
      //.selectAll(".box")
        .transition().duration(animation_speed/2)
        .style("fill", function(d){
          return coloring[d.cancer];
        });

      if (boolText){
        selection.selectAll("text")
          .transition().duration(animation_speed/2)
          .style("font-size", (boxWidth < 10) ? boxWidth : 10)
          .style("fill", "#7F8C8D");
      }
    }

    function renderBoxes(selection, gene_index){

      var keep = selection.filter(function(d, i){
        return x(i) >= (labelWidth) && x(i) <= (graphWidth);
      })
      .style("fill-opacity", 1)
      .style("stroke-opacity", 0);

      var fade = selection.filter(function(d, i){
        return (x(i) < (labelWidth) || x(i) > graphWidth - boxWidth);
      })
      .style("fill-opacity", 0.25)
      .style("stroke-opacity", 0);

      selection.attr("transform", function(d, i){
        return "translate(" + x(i) +  ")";
      })

      selection.selectAll("text")
        .style("font-size", (boxWidth < 10) ? boxWidth : 10)
        .attr("transform", "translate(" + boxWidth/2 +"," + boxHeight + "), rotate(90)");

      renderOncoOpacities(keep.selectAll(".box"), 0, 1)
      renderOncoOpacities(fade.selectAll(".box"), 0, 0.25)

      selection.selectAll(".inactivating")
        .filter(function(d){
          return d.inactivating;
        })
        .attr("width", boxWidth)
        .attr("y", function(d, i){
          return ((gene_index[d.gene] ? gene_index[d.gene]: 0) + 0.375)* (boxHeight + boxMargin) + labelHeight + boxMargin;
      });

      selection.selectAll(".box")
        .filter(function(d){ return (gene_display_list.indexOf(d.gene) > -1) })
        .attr("width", boxWidth)
        .attr("height", function(d, i){
          if (d.cna){
            return boxHeight/2;
          }
          return boxHeight;
        })
        .attr("y", function(d, i){
          var index = (gene_index[d.gene] ? gene_index[d.gene]: 0),
              CNA_addition = 0;
          if (d.cna == "del"){
            CNA_addition += boxHeight/2;
          }
          return (index)* (boxHeight + boxMargin)+ CNA_addition + labelHeight + boxMargin;
      });
    }

    function renderLabels(selection, index){

      var keep = selection
      .transition().duration(animation_speed)
      .attr("transform", function(d, i){
          return "translate(0 , " + (labelHeight +index[d]*(boxHeight + boxMargin)) +  ")"
        });

      keep.selectAll("rect")
        .attr("width", width - 2*boxMargin)
        .attr("height", boxHeight)
        .attr("transform", function(d, i){
          return "translate(" +boxMargin + ", " + boxMargin + ")"
        });

      keep.selectAll("text")
        .attr("transform", function(d, i){
          return "translate(" + (2*boxMargin) + ", " + (boxHeight) +  ")"
        });
    }

    function renderOncoOpacities(selection, minOpacity, maxOpacity){

      selection
        .filter(function(d){ 
          return (gene_display_list.indexOf(d.gene) <= -1) 
        })
        .style("fill-opacity", minOpacity)
      selection
        .filter(function(d){ 
          return (gene_display_list.indexOf(d.gene) > -1) 
        })
        .style("fill-opacity", maxOpacity)
    }

    function checkValidGene(data, list){
      var valid = false;
      list.forEach(function(d){
        for (var i = 0; i<data.length; i++){
          if (d == data[i]["gene"]){
            valid = true;
            break;
          }            
        }
      })
      return valid;
    }
  }


















  // INPUT: List of [{genes: [{gene: name}, ...]}, { ... }, ...]
  function createGeneOrder(sample){
    var sort_dict = countGenes(sample, "gene");
        genes_list = d3.keys(sort_dict);
    genes_list.sort(function(a, b){
      return sort_dict[b] - sort_dict[a]
    });

    var return_dict = {},
        counter = 0,
        leftover_array = [];
    for (var i = 0; i < genes_list.length; i++){
      if (gene_display_list.indexOf(genes_list[i]) > -1){
        return_dict[genes_list[i]] = counter;
        counter++;
      }
      else{
        leftover_array.push(genes_list[i])
      }
    }
    for (var i = 0; i < leftover_array.length; i++){
      return_dict[leftover_array[i]] = counter;
      counter++;
    }
    return return_dict
  }

  // INPUT: List of [{genes: [{gene: name}, ...]}, { ... }, ...]
  function countGenes(sample, category){
    var sort_dict = {"total": 0}
    sample.forEach(function(d){
      var gene_list = d.genes;
      gene_list.forEach(function(gene){
        if (!(gene[category] in sort_dict)){
          sort_dict[gene[category]] = 0;
        }
        sort_dict[gene[category]] += 1;
        sort_dict["total"] += 1;
      });
    });
    return sort_dict;
  }

  // INPUT: List of [{genes: [{gene: name}, ...]}, { ... }, ...]
  //        Category e.g. "cancer"
  function countCategory(sample, category){
    var sort_dict = {"total": 0}
    sample.forEach(function(d){
      if (!(d[category] in sort_dict)){
        sort_dict[d[category]] = 0;
      }
      sort_dict[d[category]] += 1;
      sort_dict["total"] += 1;
    });
    return sort_dict;
  }

  function sortViaCategory(samples, category){
    var sort_dict = countCategory(samples, category);
    var buckets = {};
    var manifest = d3.keys(sort_dict).sort(
      function(a, b){
        return sort_dict[a] - sort_dict[b]
    });
    manifest.forEach(function(bucket){
      buckets[bucket] = []
    });
    samples.forEach(function(sample){
      buckets[sample[category]].push(sample);
    });
    var return_list = [];
    for (var i = 0; i < manifest.length; i ++){
      return_list = return_list.concat(buckets[manifest[i]]);
    }
    return return_list;

  }

  // INPUT: List of [{genes: [{gene: name}, ...]}, { ... }, ...]
  //        List of [gene1, gene2, gene3 ... ]
  //        int
  function sortViaCoocurrence(sample, gene_order, depth){
    if (gene_order.length == 0 || sample.length ==0){
      return sample;
    }
    else{
      var buckets = { "present_x": [],
                      "present_nx": [],
                      "not-present": []},
          copied_gene_order = gene_order.slice(),
          cur_gene = copied_gene_order.shift();
      
      sample.forEach(function(d){
        var cur_gene_list = d.genes,
            present = false;

        for (var i = 0; i < cur_gene_list.length; i++){
          if (cur_gene_list[i]["gene"] == cur_gene){
            if (depth == cur_gene_list.length){
              buckets["present_x"].push(d);
            }
            else{
              buckets["present_nx"].push(d);
            }
            present = true;
          }
        }
        if (!present){
          buckets["not-present"].push(d);
        }
      });
      var new_depth = depth +1,
          present_x = buckets["present_x"].slice(),
          recur_present_nx = sortViaCoocurrence(buckets["present_nx"].slice(), copied_gene_order.slice(), new_depth),
          recur_not_present = sortViaCoocurrence(buckets["not-present"].slice(), copied_gene_order.slice(), depth);
      return present_x.concat(recur_present_nx, recur_not_present);
    }
  }
  function makeInteractive(dashboard, selection, accessor, category){
    selection.on("mouseover", function(d){
      dashboard.highlightCategory(accessor(d), category)
    })
    .on("mouseout", function(d){
      dashboard.unhighlightCategory()
    })
    .on("click", function(d){
      dashboard.toggleCategory(accessor(d), category)
    })
  }

  function highlightGenes(selection, accessor, gene, attribute, highlightState){
      selection.filter(function(d){
        return accessor(d) == gene;
      })
      .style(attribute, highlightState);
  }

  function unhighlightGenes(selection, accessor,  attribute, normalState, hiddenState){
    
    selection.filter(function(d){
      return gene_display_list.indexOf(accessor(d)) > -1;
    })
    .style(attribute, normalState);

    selection.filter(function(d){
      return gene_display_list.indexOf(accessor(d)) <= -1;
    })
    .style(attribute, hiddenState);
  }

});