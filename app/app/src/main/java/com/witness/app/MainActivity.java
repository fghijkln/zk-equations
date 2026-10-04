package com.witness.app;

import android.app.Activity;
import android.os.Bundle;
import android.webkit.WebSettings;
import android.webkit.WebView;
import android.webkit.WebViewClient;

import com.chaquo.python.Python;
import com.chaquo.python.android.AndroidPlatform;

/**
 * Witness: 零知识方程证明。
 * WebView 加载本地 file:// 页面，JS 经 ZkBridge 直调 Chaquopy Python 引擎。
 * witness（秘密解）只在内存里参与证明生成，从不离开本机、不进证明。
 */
public class MainActivity extends Activity {

    private WebView wv;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        if (!Python.isStarted()) {
            Python.start(new AndroidPlatform(this));
        }

        wv = findViewById(R.id.webview);
        WebSettings s = wv.getSettings();
        s.setJavaScriptEnabled(true);
        s.setDomStorageEnabled(true);
        s.setAllowFileAccess(true);
        wv.addJavascriptInterface(new ZkBridge(), "WitnessApi");
        wv.setWebViewClient(new WebViewClient());
        wv.loadUrl("file:///android_asset/www/index.html");
    }

    @Override
    public void onBackPressed() {
        if (wv != null && wv.canGoBack()) {
            wv.goBack();
        } else {
            super.onBackPressed();
        }
    }
}
