// -*- JavaScript -*-
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                                   Jiao Lin
//                      California Institute of Technology
//                       (C) 2008-2009 All Rights Reserved  
//
// {LicenseText}
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//


// requires:
//    * luban-core.js


(function(luban, $) {

  // aliases
  var ef = luban.elementFactory;
  var widgets = luban.widgets;


  // table
  //  factory
  ef.table = function(kwds, parent, docmill) {
    // 
    Date.firstDayOfWeek = 7;
    Date.format = "mm/dd/yyyy";

    // the table container
    //var thetablediv = tag('div', {id: kwds.id+'-container'});
    //parent._je.append(thetablediv);

    // the table
    var column_descriptors = createColumnDescriptors(kwds);
    //var thetable = tableFactory.createTable(kwds.id, thetablediv, column_descriptors);
    var thetable = tableFactory.createTable(kwds.id, parent._je, column_descriptors);

    // the columns that act as row identifier
    var row_identifiers = kwds.model.row_identifiers
    if (row_identifiers != null) {
      thetable.data("row-identifying-cols", row_identifiers);
    }

    // data
    var data = kwds.data, rows = [];
    for (var i in data) {
      var row = {
	'id': i,
	'data': data[i]
      };
      rows.push(row);
    }
    thetable.table_appendrows_dataonly(rows);

    // handle cell-changed event
    var oncellchanged = kwds.oncellchanged;
    var oncellchanged_handler = null;
    if (oncellchanged != null && oncellchanged!='') {
      //jscode = self.compile(oncellchanged);
      oncellchanged_handler = function(cell) {
	thetable.setRowChangedDataFromCell(cell);
	docmill.compile(oncellchanged);
	return false;
      };
    }
    
    // for each column that is editable, set that
    if (kwds.view.editable) {
      for (var i in column_descriptors) {
	var coldescriptor = column_descriptors[i];
	if (coldescriptor.editable == null || !coldescriptor.editable) continue;
	var name = coldescriptor.name;
	thetable.find("td[colname='"+name+"']").dblclick(function() {
	    $(this).enable_cell_editing(oncellchanged_handler);
	  });
      }
    }
    thetable.addClass( "tabulated" );
    var kls = kwds.Class;
    if (kls) thetable.addClass(kls);
    
    return thetable.lubanElement('table');
    return thetablediv.lubanElement('table');
  };
  //  object
  widgets.table = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.table.prototype = new widgets.base;
  //
  // extensions to jquery needed for table
  // set rowChangedData of the table from the cell
  $.fn.setRowChangedDataFromCell = function(cell) {
    var table = $(this);
    var cell = $(cell);
    var row = cell.parent();
    
    var data = {};
    var row_identifying_cols = table.data('row-identifying-cols');
    for (var i in row_identifying_cols) {
      var colname = row_identifying_cols[i];
      var value = row.get_colvalue_from_row(colname);
      data[colname] = value;
    }
    data[cell.attr('colname')] = cell.extract_data_from_cell();

    table.data('row-changed-data', data);
  };




  // helpers

  // from table description, create column descriptors
  function createColumnDescriptors(table) {

    // model
    var model = table.model;
    var measures = model; // this might change later to a more complex structure

    // view
    var view = table.view;
    var cols = view.columns;
    
    // now create column descriptors
    var column_descriptors = {};
    for (var i in cols) {
      var col = cols[i];
      var measurename = col.measure;
      var measure = measures[measurename];
      var name = measure.name;
      var descriptor = column_descriptors[name] = {};

      descriptor['name'] = name
      descriptor['text'] = col.label;
      descriptor['editable'] = col.editable && view.editable;
      descriptor['datatype'] = descriptortype2tabletype[measure.type];
    }

    return column_descriptors;
  };


  function tag(name, kwds) {

    var assignments = [];
    
    for (var key in kwds) {
      var value = kwds[key];
      assignments.push( key + '=' + '"' + value + '"' );
    }
    
    var  s = "<" + name + ' ' + assignments.join(' ') + ">" + "</"+name+">";
    
    return $(s);
  };

  // prepend 'actor.' to keys
  function prependActor(s) {
    return 'actor.'+s;
  }

  // map descriptor type to table type system
  descriptortype2tabletype = {
    'str': 'text',
    'date': 'date',
    'link': 'link'
  };

 })(luban, jQuery);


// End of file
