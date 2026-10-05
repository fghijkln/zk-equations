package com.witness.app;

import android.webkit.JavascriptInterface;

import com.chaquo.python.PyObject;
import com.chaquo.python.Python;

/**
 * JS <-> Python 桥：只暴露 prove / verify 两个方法。
 * precision: 精确度分母 k（最终精度 1/k）；精确方程可传空串。
 * 返回值一律是 JSON 字符串：{"ok":true,...} 或 {"ok":false,"error":"..."}。
 */
public class ZkBridge {

    private final PyObject api;

    public ZkBridge() {
        this.api = Python.getInstance().getModule("zk_api");
    }

    @JavascriptInterface
    public String prove(String equation, String witness, String precision) {
        try {
            PyObject r = api.callAttr("prove", equation, witness, precision);
            return r.toString();
        } catch (Exception e) {
            return "{\"ok\":false,\"error\":\"" + escape(String.valueOf(e)) + "\"}";
        }
    }

    @JavascriptInterface
    public String verify(String equation, String proofJson, String precision) {
        try {
            PyObject r = api.callAttr("verify", equation, proofJson, precision);
            return r.toString();
        } catch (Exception e) {
            return "{\"ok\":false,\"error\":\"" + escape(String.valueOf(e)) + "\"}";
        }
    }

    @JavascriptInterface
    public String solve(String equation, String decimals) {
        try {
            PyObject r = api.callAttr("solve", equation, decimals);
            return r.toString();
        } catch (Exception e) {
            return "{\"ok\":false,\"error\":\"" + escape(String.valueOf(e)) + "\"}";
        }
    }

    private static String escape(String s) {
        return s.replace("\\", "\\\\")
                .replace("\"", "\\\"")
                .replace("\n", " ")
                .replace("\r", " ");
    }
}
