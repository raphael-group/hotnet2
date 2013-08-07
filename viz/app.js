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
      highlightColor = "#F1C40F",
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

  var networkDimensions = { "width": 500, "height" :550},
      oncoprintDimensions = { "width": 500, "height" :800},
      annotationDimensions = { "width": 500, "height": 500};

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
      .on("click", function(d, i){ 
        toggleButton(DOMcontainer, canvas, data); 
      })
      .text("...");
  }

  function toggleButton(DOMcontainer, canvas, data_set){
    if (canvas.classed("unselected")){
      canvas.transition().duration(animation_speed/2)
        .style("height", function(){
          return oncoprintDimensions["height"] + 500 + "px"
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

      DOMcontainer.selectAll("li")
        .filter(function(d){ return gene_display_list.indexOf(d.gene) > -1 })
        .transition().duration(animation_speed/2)
        .style("color", textColorLightest)
      DOMcontainer.selectAll("li")
        .filter(function(d){ return gene_display_list.indexOf(d.gene) <= -1 })
        .transition().duration(animation_speed/2)
        .style("color", blockColorLight)
//      loliplot.updateGenes(updated_list);
    }

    highlightCategory = function highlightCategory(datum, category){
      network.highlightGenes(datum);
      oncoprint.highlightGenes(datum);
      DOMcontainer.selectAll("li")
        .filter(function(d){ return d.gene == datum })
        .transition().duration(animation_speed/2)
        .style("background-color", blockColorStrongest)
//      loliplot.highlightGenes(updated_list);
    }

    unhighlightCategory = function unhighlightCategory(){
      network.unhighlightGenes();
      oncoprint.unhighlightGenes();
      DOMcontainer.selectAll("li")
        .transition().duration(animation_speed/2)
        .style("background-color", blockColorMedium)
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
    var network = new generateNetwork(canvas.append("div").attr("class", "network"), data_set, this);
    var oncoprint = new generateOncoprint(canvas.append("div").attr("class", "oncoprint"), data_set, this);
    var loliplot = new generateLoliplot(canvas.append("div").attr("class", "loliplot"), data_set, this);

    canvas.attr("class", "selected");


    makeInteractive(this, DOMcontainer.selectAll("li"), function(d){ return d.gene; }, "gene");

    oncoprint.updateGenes(display_dict["gene"].slice());

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
      .attr("class", "node")
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

      hideNodes.selectAll(".node")  
        .transition().duration(animation_speed)
        .attr("r", 5);

      var keepNodes = node.filter(function(d){
          return gene_display_list.indexOf(d.gene) > -1;
        });

      keepNodes.style("fill-opacity", 1);

      keepNodes.selectAll(".node")
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
      .selectAll(".node")
//      .transition().duration(animation_speed)
      .attr("r", function(d, i){
        var tempRad = (gene_display_list.indexOf(d.gene) > -1) ? count_dict[d.gene]/count_dict["total"]*100 : 0;
        return (tempRad > 15) ? tempRad : 15
      });

      var highlight = node.filter(function(d, i){
          return d.gene == gene && gene_display_list.indexOf(gene) > -1;
        });

      highlight.selectAll(".node")
        .transition().duration(animation_speed/2)
        .style("stroke", highlightColor)
        .style("stroke-width", 5)
        .style("stroke-opacity", 1);
    }

    function unhighlightGenes(){
      
      node.selectAll(".node").style("stroke-width", 0);

      node
        .selectAll(".node")
//        .transition().duration(animation_speed)
        .attr("r", function(d, i){
          return (gene_display_list.indexOf(d.gene) > -1) ? count_dict[d.gene]/count_dict["total"]*100 : 5;

        })
        .style("stroke-opacity", 1);
    }
  }

  function generateLoliplot(DOMcontainer, data_set, dashboard){
    var width = annotationDimensions["width"],
        height = annotationDimensions["height"];

    var svg = DOMcontainer
      .append("svg")
      .attr("width", width)
      .attr("height", height);
  }

  function generateOncoprint(DOMcontainer, data_set, dashboard){

    var width = oncoprintDimensions["width"],
        height = oncoprintDimensions["height"],
        labelHeight = 0,
        labelWidth = 150
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

    var zoom = d3.behavior.zoom()
      .x(x)
      .scaleExtent([1, Math.round(20*numSamples/graphWidth)])
      .on("zoom", function(){
        renderOncoprint();
      });

    var svg = DOMcontainer
      .append("svg")
      .attr("width", width)
      .attr("height", height)
      .call(zoom);

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
    var sorted_samples = sortViaCoocurrence(samples, 
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

    data_box.selectAll(".oncosample")
      .append("text")
      .attr("text-anchor", "end")
      .text(function(d){
        var nameLen = d.name.length;
        return d.name.slice(nameLen - 4, nameLen);
      })
      .attr("fill", blockColorMedium);

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

    this.renderOncoprint = renderOncoprint;
    this.highlightGenes = highlightGenes;
    this.unhighlightGenes = unhighlightGenes;
    this.updateGenes = updateGenes;

    function renderOncoprint(){
      updateDimensions();
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
      selection.selectAll(".box")
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
      selection.selectAll(".box")
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
        return (x(i) < (labelWidth) || x(i) > (graphWidth));
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

      selection.selectAll(".box")
        .filter(function(d){ return (gene_display_list.indexOf(d.gene) > -1) })
        .attr("width", boxWidth)
        .attr("height", function(d, i){
          if (d.CNA){
            return boxHeight/2;
          }
          return boxHeight;
        })
        .attr("y", function(d, i){
          var index = (gene_index[d.gene] ? gene_index[d.gene]: 0),
              CNA_addition = 0;
          if (d.CNA == "del"){
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

});