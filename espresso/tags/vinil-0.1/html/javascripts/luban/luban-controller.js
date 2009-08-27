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

    $.get(url, allargs, callback, responsetype);

    return;
  };

  // given instructions to change user interface, execute them
  C.execUIUpdateInstructions = function(data, textStatus) {

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


  // clear error boxes in the form
  C.clearFormErrorAlerts = function(form) {
    form.find('.error').hide().text('');
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
