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


// controller object
C = luban.Controller = {
  // the following two global parameters necessary for correct operations of luban

  'url': null, // controller's url: eg http://your.web.site/main.py
  'credential': {} // credential dictionary
};


// extend jquery
(function ($, luban) {

  // aliases
  var ef = luban.elementFactory;
  var C  = luban.Controller;
  var widgets = luban.widgets;

  // *** "widgets" ***
  //
  // credential
  ef.credential = function (kwds) {
    C.credential = {username: kwds.username, ticket: kwds.ticket};
    // a invisible div
    var div = $('<div></div>').hide();
    return div.lubanElement('credential');
  };
  widgets.credential = function(elem) {
    this.super1 = widgets.base;
    this.super1(elem);
  };
  widgets.credential.prototype = new widgets.base;




  // method to submit form (this) to the given actor and routine
  // actor and routine are specified in kwds
  // kwds:
  //   actor: name of actor
  //   routine: name of routine
  //   callback: callback function when response from the server is obtained
  //   responsetype: type of response obtained from server. default: json
  //   data: additional data to send to the controller
  $.fn.submitTo = function(kwds) {
    var actor = kwds.actor;
    var routine = kwds.routine;
    if (!routine) routine = 'default';
    var callback = kwds.callback;
    var controller = C.url;

    var responsetype = kwds.responsetype;
    if (responsetype == null)
      responsetype ='json';

    var formdatastr = $(this).serialize();

    var data = kwds.data, datastr;
    if (data == null) datastr = '';
    else datastr = argStr(data);

    var sentrystr = argStr(C.getCredentialArgs());
    
    var args = {'actor': actor,
		'routine': routine};
    var allargsstr = [argStr(args), sentrystr, datastr, formdatastr].join('&');

    // jquery ajax post
    $.post(controller, allargsstr, callback, responsetype);
    
    return $(this);
  };

  
  // replace the content of "this" widget
  // data: a dict
  //   html: new html content
  //   includes: paths of js scripts to include
  //   script: js script to run
  $.fn.replaceContent = function(data) {
    // clear
    $(this).empty();

    C.execUIUpdateInstructions(data);
  };


  // helper function
  function argStr(args) {
    var assignments = [];
    for (var k in args) {
      v = args[k];
      assignment = k+'='+v;
      assignments.push(assignment);
    }
    return assignments.join('&');
  }
  
  // prepend 'actor.' to keys
  function prependActorStr(args) {
    var d = {};
    for (var k in args) {
      var k1 = 'actor.'+k;
      d[k1] = args[k];
    }
    return d;
  }


  // controller methods

  // call controller
  // kwds
  //   actor: the name of the actor
  //   routine: the name of the routine
  //   callback: the call back function when the response of the server is received
  //   responsetype: the expected response type. default: json
  //   data: the additional data to send to the server
  C.call = function (kwds) {
    var actor = kwds.actor;
    var routine = kwds.routine;
    if (!routine) routine='default';
    var callback = kwds.callback;
    var url = C.url;

    var responsetype = kwds.responsetype;
    if (responsetype==null)
      responsetype = 'json';

    // call
    var args = {'actor': actor,
		'routine': routine};

    // data
    var data = kwds.data;
    if (data == null) data = {};

    // credential
    var credArgs = C.getCredentialArgs();

    // all
    var allargs = $.extend({}, args, data, credArgs);

    C.runWithLoadingAlert(function(callback) {
	$.get(url, allargs, callback, responsetype);
      }, 
      callback);
    
    return;
  };

  // run a function (which has a callback function when the function finishes)
  // with "loading ..." alert shown on the window
  C.runWithLoadingAlert = function (func, callback) {
    var callback1 = function () {
      // shut down the loading alert
      C.notifyLoadingDivToEnd(func);

      // start callback function
      callback.apply({}, arguments);
    };
    // start the loading alert
    C.notifyLoadingDivToStart(func);
    // call the function with new callback
    func(callback1);
    return;
  };

  C.notifyLoadingDivToStart = function (func) {
    var loadingdiv = C.getLoadingDiv();
    var f = function () {
      if (loadingdiv.data('running')) {
	return;
      }
	  
      // updating the 'loading' message periodically
      var intervals = loadingdiv.data('intervals');
      intervals[func] = setInterval("$('#luban-----loadingdiv').trigger('update-loading-alert')", 1000);

      var wh=$(window).height(), ww=$(window).width();
      loadingdiv.css('left', ww/2-20);
      loadingdiv.css('top', 5); //wh/2);
      loadingdiv.show();      
    }

    setTimeout(f, 1000);
  };

  C.notifyLoadingDivToEnd = function (func) {
    var loadingdiv = C.getLoadingDiv();

    var f = function () {
      var intervals = loadingdiv.data('intervals');
      var interval = intervals[func];
      clearInterval(interval);
      delete intervals[func];
      
      loadingdiv.hide();
      loadingdiv.data('running', 0);
    };
    setTimeout(f, 1200); // this timeout must be longer than the timeout in notifyLoadingDivToStart, otherwise some interval may not be killed correctly

  };

  C.getLoadingDiv = function () {
    var id = 'luban-----loadingdiv';
    var div = $('#'+id);

    // this should not happend
    if (div.length>1) throw "?";

    // if not found, creat it
    if (div.length==0) {

      // create the div
      div = $('<div id="'+id+'"/>');
      div.text('loading...');
      div.hide();
      
      // attach to page
      $('body').append(div);

      //
      div.data('intervals', {});
      div.data('running', 0);

      // events
      div.bind('update-loading-alert', function() {
	  var $this = $(this);
	  $this.data('running', 1);
	  var ndots = $this.data('ndots');
	  if (ndots == null) ndots = 0;
	  else ndots += 1;
	  if (ndots > 3) ndots=0;
	  $this.data('ndots', ndots);
	  $this.text('loading'+Array(ndots+1).join('.'));
	});
    }
    return div;
  };

  C.getErrorReportDiv = function () {
    var id = 'luban-----errorreport';
    var div = $('#'+id);

    // this should not happend
    if (div.length>1) throw "?";

    // if not found, creat it
    if (div.length==0) {

      // create the div
      div = $('<div id="'+id+'"/>');
      
      // create a title section
      var titlediv = $('<div class="title-container"/>');
      div.append(titlediv);
      // create an interior container
      var interior_div = $('<div class="interior-container"/>');
      div.append(interior_div);

      // create a ok button
      var okdiv = $("<div class='ok-container'/>"); div.append(okdiv);
      var ok = $('<input type="submit" value="OK"/>'); okdiv.append(ok);
      ok.click(function () {
	  var div = $('#luban-----errorreport');
	  div.hide();
	  div.find('.interior-container').empty();
	});

      // attach to page
      $('body').append(div);
    }

    return div;
  };

  // given instructions to change user interface, execute them
  C.execUIUpdateInstructions = function(data, textStatus) {
    
    var exception = data.exception;
    if (exception) {
      var div = C.getErrorReportDiv();
      var width = $(window).width();
      var height = $(window).height();
      div.css('left', 50);
      div.css('top', 20);
      div.width(width-100);
      //div.css('max-height', height-60);
      div.height(height-60);
      var title = div.children('.title-container');
      title.html('<h3>Exception raised</h3>');
      var interior = div.children('.interior-container');
      interior.html('<pre>'+exception+'</pre>');
      div.show();
      return;
    }
    
    // decompose data
    var html = data.html;
    var includes = data.includes;
    var script = data.script;
    
    // don't know how html can be useful at this moment
    // html;
    
    // include scripts
    var commands = [];
    for (var index in includes) {
      var include = includes[index];
      var f = function (callback) {
	$.getScript(include, callback);
      };
      commands.push(f);
    };
    
    // run script
    commands.push( function () {eval(script);} );
    
    runCmds(commands);
  };


  // run a sequence of commands
  // each command is a javascript function that has the signature
  // function command(callback)
  // callback will be called when command is done.
  function runCmds(cmds) {
    if (cmds.length==1) return cmds[0]();
    return cmds[0]( runCmds(cmds.slice(1)) );
  }


  // load from server and execute commands in the response
  // kwds: a dict
  //   actor: the name of the actor
  //   routine: the name of the routine
  //   data: a dictionary of additional parameters to send to the server
  C.load = function(kwds) {
    var data = kwds.data;
    kwds.data = prependActorStr(data);
    var callback = C.execUIUpdateInstructions;
    kwds.callback = callback;
    C.call(kwds);
  };


  // notify the server the event happened to an element,
  // and execute commands in the response
  // element: the element where event happened
  // event: the name of the event
  // kwds: a dict
  //   actor: the name of the actor
  //   routine: the name of the routine
  //   data: a dictionary of additional parameters to send to the server
  C.notify = function(element, event, kwds) {
    var evtdata = element.getEventData(event);
    var data = kwds.data;
    var tmp = $.extend({}, data, evtdata);
    tmp = prependActorStr(tmp);
    kwds.data = tmp;
    
    var callback = C.execUIUpdateInstructions;
    kwds.callback = callback;
    C.call(kwds);
  };


  // submit form to server and get response and execute commands in the response
  // form: the form to submit
  // kwds: a dict
  //   actor: the name of the actor
  //   routine: the name of the routine
  //   data: a dictionary of additional parameters to send to the server
  C.submit = function(form, kwds) {
    var callback = C.execUIUpdateInstructions;
    kwds.callback = callback;
    kwds.data = prependActorStr(kwds.data);
    form.submitTo(kwds);
  };


  // create the credentail data to be send to the server
  C.getCredentialArgs = function() {
    var credential = C.credential;
    if (credential == null)
      return {};

    var ret = {};
    var username = credential.username;
    if (username != null) ret['sentry.username'] = username;
    
    var ticket = credential.ticket;
    if (ticket != null) ret['sentry.ticket'] = ticket;

    return ret;
  };
  
 })(jQuery, luban);


// End of file
