/* Javascript for MultiQC Default Template */

var brewer_scales = ['￼YlOrRd', 'YlOrBr', 'YlGnBu', 'YlGn', 'Reds', 'RdPu',
  'Purples', 'PuRd', 'PuBuGn', 'PuBu', 'OrRd', 'Oranges', 'Greys', 'Greens',
  'GnBu', 'BuPu', 'BuGn', 'Blues', '￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼Set3', 'Set2', 'Set1', '￼￼￼Pastel2', 'Pastel1',
  'Paired', 'Dark2', 'Accent', '￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼Spectral', 'RdYlGn', 'RdYlBu', 'RdGy', 'RdBu',
  'PuOr', 'PRGn', 'PiYG', 'BrBG'];

// Execute when page load has finished
$(function () {

  // Enable the bootstrap tooltip hovers
  // http://getbootstrap.com/javascript/#tooltips
  $('[data-toggle="tooltip"]').tooltip();

  // Hide suspected duplicates in general stats table
  hide_general_stats_duplicates();

  $('#genstat_table_showhide_dups').click(function(e){
    e.preventDefault();
    if($(this).attr('data-ishidden') === 'true'){
      $('.genstats-duplicate').show();
      $(this).text('Click here to hide suspected duplicate samples.').attr('data-ishidden', 'false');
    } else {
      $('.genstats-duplicate').hide();
      $(this).text('Click here to show suspected duplicate samples.').attr('data-ishidden', 'true');
    }
  });


  // Colour code table cells using chroma.js
  $('table').each(function(){
    var table = $(this);
    table.find('thead th').each(function(index){
      if($(this).hasClass('chroma-col')){

        // Get the colour scheme if set
        var colscheme_rev = false;
        var colscheme = $(this).data('chroma-scale');
        if(colscheme.substr(colscheme.length - 4) == '-rev'){
          colscheme_rev = true;
          colscheme = colscheme.substr(0, colscheme.length - 4);
        }
        if(colscheme === undefined || brewer_scales.indexOf(colscheme) == -1){
          colscheme = 'GnBu';
        }

        // Collect the data series
        var data = [];
        table.find('tr td:nth-of-type('+(index)+')').each(function(){
          var d = parseFloat($(this).text());
          data.push(d);
        });
        if(data.length == 0){ return true; } // Skip to next loop

        // Get the max and min values if not set with data attributes
        var maxval = $(this).data('chroma-max');
        var minval = $(this).data('chroma-min');
        if(maxval === undefined || minval === undefined){
          $.each(data, function(k, v){
            if(v > maxval || maxval == undefined){ maxval = v; }
            if(v < minval || minval == undefined){ minval = v; }
          });
        }
        if(isNaN(minval) || isNaN(maxval)){
          console.log('Could not calculate max or min value for '+$(this).text()+': ['+[minval, maxval]+']')
          return true; // Skip to next loop
        }

        // Go through table cells again, adding colour
        var i = 0;
        var scale = chroma.scale(colscheme).domain([minval, maxval]);
        if(colscheme_rev){
          scale = chroma.scale(colscheme).domain([maxval, minval]);
        }
        table.find('tr td:nth-of-type('+(index)+')').each(function(){
          var col = scale(parseFloat($(this).text())).css();
          $(this).css('background-color', col);
          // Change the colour of the text if we can get a better contrast ratio
          if(chroma.contrast('#EEEEEE', col) > chroma.contrast($(this).css('color'), col)){
            $(this).css('color', '#EEEEEE').css('font-weight', '200');
          }
        });
      }
    });
  });


})


////////////////////////////////////////////////
// HighCharts Plotting Functions
////////////////////////////////////////////////

// Basic Line Graph
function plot_xy_line_graph(div, data, config){
  if(config['tt_label'] === undefined){ config['tt_label'] = '{point.x}'; }
  if(config['use_legend'] === undefined){ config['use_legend'] = true; }
  if(config['click_func'] === undefined){ config['click_func'] = function(){}; }
  $(div).highcharts({
    chart: {
      type: 'line',
      zoomType: 'x'
    },
    title: {
      text: config['title'],
      x: -20 //center
    },
    xAxis: {
      title: {
        text: config['xlab']
      },
      max: config['xmax'],
      min: config['xmin']
    },
    yAxis: {
      title: {
        text: config['ylab']
      },
      max: config['ymax'],
      min: config['ymin'],
      plotLines: [{
        value: 0,
        width: 1,
        color: '#808080'
      }]
    },
    plotOptions: {
      series: {
        marker: { enabled: false },
        cursor: 'pointer',
        point: {
          events: {
            click: config['click_func']
          }
        }
      }
    },
    legend: {
      enabled: config['use_legend'],
      layout: 'vertical',
      align: 'right',
      verticalAlign: 'middle',
      borderWidth: 0
    },
    credits: {
			enabled: false
		},
    tooltip: {
      headerFormat: '',
			pointFormat: '<span style="color:{series.color}; text-decoration:underline;">{series.name}</span><br>'+config['tt_label']+': {point.y:.2f}',
			useHTML: true
    },
    series: data
  });
}

// Hide suspected duplicates in general stats table
function hide_general_stats_duplicates(){
  var gen_stats_snames = [];
  var hidden_rows = 0;
  $("#general_stats_table tbody tr").each(function(){
    var sn = $(this).find('th').text();
    var matched = 0;
    $.each(gen_stats_snames, function(k, v){
      if(sn.indexOf(v) > -1){
        matched += 1;
      }
    });
    if(matched > 0){
      $(this).addClass('genstats-duplicate').hide();
      hidden_rows += 1;
    } else {
      gen_stats_snames.push(sn)
    }
  });
  if(hidden_rows > 0){
    $('#genstat_table_showhide_dups').text('Click here to show suspected duplicate samples.').show();
  }
}



//////////////////
// Generic helper functions
// http://stackoverflow.com/questions/6735470/get-pixel-color-from-canvas-on-mouseover
function findPos(obj) {
  var curleft = 0, curtop = 0;
  if (obj.offsetParent) {
    do {
      curleft += obj.offsetLeft;
      curtop += obj.offsetTop;
    } while (obj = obj.offsetParent);
    return { x: curleft, y: curtop };
  }
  return undefined;
}
