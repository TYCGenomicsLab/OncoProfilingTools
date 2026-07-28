(function(){
  function onReady(callback){if(document.readyState==='loading'){document.addEventListener('DOMContentLoaded',callback);}else{callback();}}
  onReady(function(){
    if(!window.Shiny)return;
    Shiny.addCustomMessageHandler('module-compatibility',function(message){
      var card=document.getElementById('module-card-'+message.key);
      var badge=document.getElementById('compat-'+message.key);
      var input=document.getElementById('module_'+message.key);
      if(!card||!badge||!input)return;
      var compatible=Boolean(message.compatible&&message.ready);
      card.classList.toggle('module-incompatible',message.ready&&!compatible);
      badge.classList.toggle('compatible',compatible);
      badge.textContent=!message.ready?'Awaiting data':compatible?'Compatible':'Input required';
      input.disabled=!compatible;
      if(message.ready){input.checked=compatible;input.dispatchEvent(new Event('change',{bubbles:true}));}
    });
    Shiny.addCustomMessageHandler('analysis-progress',function(message){
      var row=document.getElementById('progress-'+message.key);
      if(!row)return;
      row.className='progress-row progress-'+message.status;
      var icon=row.querySelector('.progress-indicator');
      var copy=row.querySelector('.progress-message');
      if(icon)icon.textContent=message.status==='completed'?'✓':message.status==='error'?'!':message.status==='running'?'↻':message.status==='skipped'?'–':'·';
      if(copy)copy.textContent=message.message||message.status;
    });
    Shiny.addCustomMessageHandler('reset-analysis-ui',function(){
      document.querySelectorAll('.progress-row').forEach(function(row){
        row.className='progress-row progress-idle';
        var icon=row.querySelector('.progress-indicator');
        var copy=row.querySelector('.progress-message');
        if(icon)icon.textContent='·';
        if(copy)copy.textContent='Waiting';
      });
    });
  });
})();
